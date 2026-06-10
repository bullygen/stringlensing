from __future__ import annotations

import csv
import math
import os
import re
import shlex
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

os.environ.setdefault("MPLCONFIGDIR", str(Path("/tmp") / "matplotlib-eda-photometry"))

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from photutils.detection import DAOStarFinder
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QApplication,
    QCheckBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QSplitter,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)
from scipy.ndimage import median_filter
from scipy.optimize import curve_fit


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DEFAULT_INPUT_DIR = PROJECT_ROOT / "targets_photometry"
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "SL_photometry_real"

GAIN_KEYS = ("GAIN", "EGAIN", "CCDGAIN", "ATODGAIN", "ATODGNA")
MANIFEST_COLUMNS = (
    "sample_id",
    "image_path",
    "truth_json_path",
    "env_nx",
    "env_ny",
    "env_pixel_scale",
    "env_noise_sigma",
    "env_psf_sigma",
)


@dataclass
class ImageSlice:
    file_path: Path
    hdu_index: int
    hdu_name: str
    slice_index: tuple[int, ...]
    label: str
    data: np.ndarray
    header: fits.Header


@dataclass
class LoadedFile:
    file_path: Path
    status: str
    slices: list[ImageSlice] = field(default_factory=list)
    message: str = ""


@dataclass
class BackgroundResult:
    image_electrons: np.ndarray
    background_mean_adu: float
    background_sigma_adu: float
    background_mean_electrons: float
    background_sigma_electrons: float
    gain: float
    gain_source: str
    masked_fraction: float
    mask: np.ndarray


@dataclass
class PsfResult:
    n_sources: int
    n_fitted: int
    sigma_pix: float
    fwhm_arcsec: float


def safe_sample_id(path: Path, label: str) -> str:
    raw = f"{path.stem}_{label}"
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", raw)
    return cleaned.strip("_") or "sample"


def finite_image(data: np.ndarray) -> np.ndarray:
    image = np.asarray(data, dtype=np.float64)
    if image.ndim != 2:
        raise ValueError(f"Expected a 2D image, got shape {image.shape}")
    return image


def display_limits(image: np.ndarray) -> tuple[float, float]:
    finite = image[np.isfinite(image)]
    if finite.size == 0:
        return 0.0, 1.0
    vmin, vmax = np.nanpercentile(finite, [1, 99])
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
        vmin = float(np.nanmin(finite))
        vmax = float(np.nanmax(finite))
    if vmin == vmax:
        vmax = vmin + 1.0
    return float(vmin), float(vmax)


def iter_displayable_slices(file_path: Path) -> LoadedFile:
    slices: list[ImageSlice] = []
    try:
        with fits.open(file_path, memmap=False, output_verify="silentfix") as hdul:
            hdul.verify("silentfix")
            primary_header = hdul[0].header.copy() if len(hdul) else fits.Header()
            for hdu_index, hdu in enumerate(hdul):
                data = hdu.data
                if data is None:
                    continue
                header = primary_header.copy()
                header.update(hdu.header)
                array = np.asarray(data)
                if array.ndim < 2:
                    continue
                if array.ndim == 2:
                    if _is_numeric_image(array):
                        label = f"HDU {hdu_index} {hdu.name}"
                        slices.append(
                            ImageSlice(
                                file_path=file_path,
                                hdu_index=hdu_index,
                                hdu_name=hdu.name,
                                slice_index=(),
                                label=label,
                                data=finite_image(array),
                                header=header,
                            )
                        )
                    continue

                leading_shape = array.shape[:-2]
                for leading_index in np.ndindex(leading_shape):
                    plane = array[leading_index]
                    if not _is_numeric_image(plane):
                        continue
                    label = f"HDU {hdu_index} {hdu.name} slice {leading_index}"
                    slices.append(
                        ImageSlice(
                            file_path=file_path,
                            hdu_index=hdu_index,
                            hdu_name=hdu.name,
                            slice_index=tuple(int(v) for v in leading_index),
                            label=label,
                            data=finite_image(plane),
                            header=header,
                        )
                    )
    except Exception as exc:
        return LoadedFile(file_path=file_path, status="invalid", message=str(exc))

    if not slices:
        return LoadedFile(file_path=file_path, status="no-image", message="No displayable 2D image slices found.")
    return LoadedFile(file_path=file_path, status="ok", slices=slices)


def _is_numeric_image(array: np.ndarray) -> bool:
    return np.issubdtype(array.dtype, np.number) and array.ndim == 2 and min(array.shape) > 1


def find_gain(header: fits.Header, image_adu: np.ndarray, background_mask: np.ndarray) -> tuple[float, str]:
    for key in GAIN_KEYS:
        value = header.get(key)
        try:
            gain = float(value)
        except (TypeError, ValueError):
            continue
        if np.isfinite(gain) and gain > 0:
            return gain, f"header {key}"

    background = image_adu[background_mask & np.isfinite(image_adu)]
    if background.size < 16:
        return 1.0, "fallback unity; too few background pixels"

    mean_adu = float(np.nanmean(background))
    var_adu = float(np.nanvar(background, ddof=1))
    if np.isfinite(mean_adu) and np.isfinite(var_adu) and mean_adu > 0 and var_adu > 0:
        return mean_adu / var_adu, "estimated mean/variance"
    return 1.0, "fallback unity; mean/variance unusable"


def estimate_background_and_gain(image_adu: np.ndarray, header: fits.Header) -> BackgroundResult:
    image = finite_image(image_adu)
    finite = np.isfinite(image)
    if finite.sum() < 16:
        raise ValueError("Image has too few finite pixels for background estimation.")

    work = np.array(image, dtype=np.float64, copy=True)
    fill_value = float(np.nanmedian(work[finite]))
    work[~finite] = fill_value
    background_mask = finite.copy()

    for _ in range(5):
        filtered = median_filter(work, size=7, mode="nearest")
        residual = image - filtered
        _, med_resid, std_resid = sigma_clipped_stats(residual[background_mask], sigma=3.0, maxiters=5)
        if not np.isfinite(std_resid) or std_resid <= 0:
            break
        background_mask = finite & (residual < med_resid + 3.0 * std_resid) & (residual > med_resid - 5.0 * std_resid)

    if background_mask.sum() < 16:
        background_mask = finite

    _, background_mean_adu, background_sigma_adu = sigma_clipped_stats(image[background_mask], sigma=3.0, maxiters=5)
    gain, gain_source = find_gain(header, image, background_mask)
    image_electrons = image * gain
    masked_fraction = 1.0 - float(background_mask.sum()) / float(finite.sum())

    return BackgroundResult(
        image_electrons=image_electrons,
        background_mean_adu=float(background_mean_adu),
        background_sigma_adu=float(background_sigma_adu),
        background_mean_electrons=float(background_mean_adu * gain),
        background_sigma_electrons=float(background_sigma_adu * gain),
        gain=float(gain),
        gain_source=gain_source,
        masked_fraction=masked_fraction,
        mask=background_mask,
    )


def gaussian_2d(coords: tuple[np.ndarray, np.ndarray], amp: float, x0: float, y0: float, sx: float, sy: float, offset: float) -> np.ndarray:
    x, y = coords
    sx = max(float(sx), 1e-6)
    sy = max(float(sy), 1e-6)
    model = offset + amp * np.exp(-0.5 * (((x - x0) / sx) ** 2 + ((y - y0) / sy) ** 2))
    return model.ravel()


def estimate_psf(image: np.ndarray, pixel_scale: float) -> PsfResult:
    data = finite_image(image)
    mean, median, std = sigma_clipped_stats(data, sigma=3.0, maxiters=5)
    if not np.isfinite(std) or std <= 0:
        raise ValueError("Cannot estimate PSF because background standard deviation is invalid.")

    finder = DAOStarFinder(fwhm=4.0, threshold=5.0 * std)
    sources = finder(data - median)
    if sources is None or len(sources) == 0:
        raise ValueError("No stars detected above the current DAOStarFinder threshold.")

    sigmas: list[float] = []
    half_size = 7
    max_sources = min(len(sources), 40)
    x_col = "xcentroid" if "xcentroid" in sources.colnames else "x_centroid"
    y_col = "ycentroid" if "ycentroid" in sources.colnames else "y_centroid"
    for row in sources[:max_sources]:
        x = int(round(float(row[x_col])))
        y = int(round(float(row[y_col])))
        y0 = max(0, y - half_size)
        y1 = min(data.shape[0], y + half_size + 1)
        x0 = max(0, x - half_size)
        x1 = min(data.shape[1], x + half_size + 1)
        stamp = data[y0:y1, x0:x1]
        if stamp.shape[0] < 7 or stamp.shape[1] < 7:
            continue
        local = stamp - np.nanmedian(stamp)
        if not np.isfinite(local).all() or np.nanmax(local) <= 0:
            continue

        yy, xx = np.indices(stamp.shape)
        p0 = [
            float(np.nanmax(local)),
            float(x - x0),
            float(y - y0),
            2.0,
            2.0,
            float(np.nanmedian(stamp)),
        ]
        lower = [0.0, 0.0, 0.0, 0.5, 0.5, float(np.nanmin(stamp))]
        upper = [
            float(max(np.nanmax(stamp) - np.nanmin(stamp), 1.0) * 4.0),
            float(stamp.shape[1] - 1),
            float(stamp.shape[0] - 1),
            10.0,
            10.0,
            float(np.nanmax(stamp)),
        ]

        try:
            popt, _ = curve_fit(
                gaussian_2d,
                (xx, yy),
                stamp.ravel(),
                p0=p0,
                bounds=(lower, upper),
                maxfev=3000,
            )
        except Exception:
            continue

        sigma = math.sqrt(abs(float(popt[3]) * float(popt[4])))
        if np.isfinite(sigma) and 0.5 <= sigma <= 10.0:
            sigmas.append(sigma)

    if not sigmas:
        raise ValueError(f"Detected {len(sources)} sources, but no Gaussian star fits were usable.")

    sigma_pix = float(np.nanmedian(sigmas))
    fwhm_arcsec = float(2.355 * sigma_pix * pixel_scale)
    return PsfResult(n_sources=int(len(sources)), n_fitted=len(sigmas), sigma_pix=sigma_pix, fwhm_arcsec=fwhm_arcsec)


class ImageCanvas(QWidget):
    def __init__(self) -> None:
        super().__init__()
        self.figure = Figure(figsize=(6, 5))
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.ax = self.figure.add_subplot(111)
        self.image_artist = None
        self.crop_patch: Optional[Rectangle] = None
        self.click_callback = None

        layout = QVBoxLayout(self)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.canvas.mpl_connect("button_press_event", self._on_click)

    def show_image(self, image: np.ndarray, title: str) -> None:
        self.ax.clear()
        vmin, vmax = display_limits(image)
        self.image_artist = self.ax.imshow(image, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)
        self.ax.set_title(title)
        self.ax.set_xlabel("x pixel")
        self.ax.set_ylabel("y pixel")
        self.crop_patch = None
        self.figure.tight_layout()
        self.canvas.draw_idle()

    def set_crop_rect(self, cx: int, cy: int, nx: int, ny: int) -> None:
        if self.crop_patch is not None:
            self.crop_patch.remove()
        x0 = cx - nx / 2.0
        y0 = cy - ny / 2.0
        self.crop_patch = Rectangle((x0, y0), nx, ny, fill=False, edgecolor="cyan", linewidth=1.6)
        self.ax.add_patch(self.crop_patch)
        self.canvas.draw_idle()

    def _on_click(self, event: Any) -> None:
        if event.inaxes != self.ax or event.xdata is None or event.ydata is None:
            return
        if self.click_callback is not None:
            self.click_callback(int(round(event.xdata)), int(round(event.ydata)))


class PhotometryEdaWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Photometry FITS EDA")
        self.resize(1450, 900)

        self.input_dir = DEFAULT_INPUT_DIR
        self.output_dir = DEFAULT_OUTPUT_DIR
        self.loaded_files: dict[Path, LoadedFile] = {}
        self.current_slice: Optional[ImageSlice] = None
        self.current_image_adu: Optional[np.ndarray] = None
        self.current_image_electrons: Optional[np.ndarray] = None
        self.background_result: Optional[BackgroundResult] = None
        self.psf_result: Optional[PsfResult] = None
        self.locked_crop_shape: Optional[tuple[int, int]] = self._read_locked_crop_shape()

        self.file_list = QListWidget()
        self.slice_list = QListWidget()
        self.canvas = ImageCanvas()
        self.canvas.click_callback = self._set_center_from_click
        self.status_box = QTextEdit()
        self.status_box.setReadOnly(True)
        self.status_box.setMinimumHeight(140)

        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(self._build_left_panel())
        splitter.addWidget(self.canvas)
        splitter.addWidget(self._build_right_panel())
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 4)
        splitter.setStretchFactor(2, 2)
        layout.addWidget(splitter)
        layout.addWidget(self.status_box)

        self._scan_default_dir()
        self._set_editing_enabled(False)
        self._log("Ready. Select FITS files and press Open selected FITS.")

    def _build_left_panel(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)

        self.input_label = QLabel(str(self.input_dir))
        self.input_label.setWordWrap(True)
        choose_dir_btn = QPushButton("Choose directory")
        choose_dir_btn.clicked.connect(self._choose_directory)
        choose_files_btn = QPushButton("Choose FITS files")
        choose_files_btn.clicked.connect(self._choose_files)
        scan_btn = QPushButton("Scan default directory")
        scan_btn.clicked.connect(self._scan_default_dir)
        open_btn = QPushButton("Open selected FITS")
        open_btn.clicked.connect(self._open_selected_files)
        discard_btn = QPushButton("Discard selected file")
        discard_btn.clicked.connect(self._discard_selected_file)

        layout.addWidget(QLabel("FITS files"))
        layout.addWidget(self.input_label)
        layout.addWidget(choose_dir_btn)
        layout.addWidget(choose_files_btn)
        layout.addWidget(scan_btn)
        layout.addWidget(open_btn)
        layout.addWidget(discard_btn)
        layout.addWidget(self.file_list, stretch=1)
        layout.addWidget(QLabel("Displayable slices"))
        layout.addWidget(self.slice_list, stretch=1)

        self.slice_list.itemSelectionChanged.connect(self._select_slice)
        return panel

    def _build_right_panel(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)

        view_group = QGroupBox("View")
        view_layout = QVBoxLayout(view_group)
        start_btn = QPushButton("Start editing")
        start_btn.clicked.connect(self._start_editing)
        view_layout.addWidget(start_btn)

        bg_group = QGroupBox("Background and gain")
        bg_layout = QVBoxLayout(bg_group)
        self.bg_btn = QPushButton("Compute background/gain")
        self.bg_btn.clicked.connect(self._compute_background)
        self.bg_label = QLabel("No background estimate yet.")
        self.bg_label.setWordWrap(True)
        bg_layout.addWidget(self.bg_btn)
        bg_layout.addWidget(self.bg_label)

        psf_group = QGroupBox("PSF")
        psf_layout = QVBoxLayout(psf_group)
        self.psf_btn = QPushButton("Estimate PSF")
        self.psf_btn.clicked.connect(self._compute_psf)
        self.psf_label = QLabel("No PSF estimate yet.")
        self.psf_label.setWordWrap(True)
        psf_layout.addWidget(self.psf_btn)
        psf_layout.addWidget(self.psf_label)

        crop_group = QGroupBox("Crop and export")
        crop_layout = QFormLayout(crop_group)
        self.pixel_scale_spin = QDoubleSpinBox()
        self.pixel_scale_spin.setDecimals(5)
        self.pixel_scale_spin.setRange(1e-5, 100.0)
        self.pixel_scale_spin.setValue(0.10)
        self.pixel_scale_spin.setSingleStep(0.01)
        self.center_x_spin = QSpinBox()
        self.center_y_spin = QSpinBox()
        self.center_x_spin.setRange(0, 1000000)
        self.center_y_spin.setRange(0, 1000000)
        self.crop_nx_spin = QSpinBox()
        self.crop_ny_spin = QSpinBox()
        self.crop_nx_spin.setRange(1, 1000000)
        self.crop_ny_spin.setRange(1, 1000000)
        self.crop_nx_spin.setValue(128)
        self.crop_ny_spin.setValue(128)
        self.output_edit = QLineEdit(str(self.output_dir))
        output_btn = QPushButton("Choose output folder")
        output_btn.clicked.connect(self._choose_output_dir)
        preview_btn = QPushButton("Preview crop")
        preview_btn.clicked.connect(self._preview_crop)
        self.save_btn = QPushButton("Save crop for fit pipeline")
        self.save_btn.clicked.connect(self._save_crop)
        self.allow_shape_reset = QCheckBox("Allow new project shape")

        crop_layout.addRow("Pixel scale arcsec/pix", self.pixel_scale_spin)
        crop_layout.addRow("Center x", self.center_x_spin)
        crop_layout.addRow("Center y", self.center_y_spin)
        crop_layout.addRow("Crop nx", self.crop_nx_spin)
        crop_layout.addRow("Crop ny", self.crop_ny_spin)
        crop_layout.addRow("Output folder", self.output_edit)
        crop_layout.addRow(output_btn)
        crop_layout.addRow(preview_btn)
        crop_layout.addRow(self.save_btn)
        crop_layout.addRow(self.allow_shape_reset)

        instruction_group = QGroupBox("Fit command")
        instruction_layout = QVBoxLayout(instruction_group)
        self.fit_command_label = QLabel()
        self.fit_command_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.fit_command_label.setWordWrap(True)
        instruction_layout.addWidget(self.fit_command_label)
        self._update_fit_command()

        layout.addWidget(view_group)
        layout.addWidget(bg_group)
        layout.addWidget(psf_group)
        layout.addWidget(crop_group)
        layout.addWidget(instruction_group)
        layout.addStretch(1)

        for widget in (self.crop_nx_spin, self.crop_ny_spin, self.pixel_scale_spin):
            widget.valueChanged.connect(lambda _value: self._update_fit_command())
        self.output_edit.textChanged.connect(lambda _text: self._update_fit_command())
        return panel

    def _set_editing_enabled(self, enabled: bool) -> None:
        for widget in (
            self.bg_btn,
            self.psf_btn,
            self.center_x_spin,
            self.center_y_spin,
            self.crop_nx_spin,
            self.crop_ny_spin,
            self.pixel_scale_spin,
            self.save_btn,
        ):
            widget.setEnabled(enabled)

    def _scan_default_dir(self) -> None:
        self.input_dir = DEFAULT_INPUT_DIR
        self.input_label.setText(str(self.input_dir))
        files = sorted(self.input_dir.rglob("*.fits")) if self.input_dir.exists() else []
        self._populate_file_list(files)
        self._log(f"Scanned {self.input_dir}: found {len(files)} FITS files.")

    def _choose_directory(self) -> None:
        directory = QFileDialog.getExistingDirectory(self, "Choose FITS directory", str(self.input_dir))
        if not directory:
            return
        self.input_dir = Path(directory)
        self.input_label.setText(str(self.input_dir))
        files = sorted(self.input_dir.rglob("*.fits"))
        self._populate_file_list(files)
        self._log(f"Scanned {self.input_dir}: found {len(files)} FITS files.")

    def _choose_files(self) -> None:
        files, _ = QFileDialog.getOpenFileNames(self, "Choose FITS files", str(self.input_dir), "FITS files (*.fits *.fit *.fts);;All files (*)")
        if files:
            self._populate_file_list([Path(path) for path in files])
            self._log(f"Selected {len(files)} FITS files from file picker.")

    def _populate_file_list(self, files: list[Path]) -> None:
        self.file_list.clear()
        self.slice_list.clear()
        self.loaded_files.clear()
        for path in files:
            item = QListWidgetItem(path.name)
            item.setToolTip(str(path))
            item.setData(Qt.ItemDataRole.UserRole, path)
            self.file_list.addItem(item)

    def _selected_file_paths(self) -> list[Path]:
        items = self.file_list.selectedItems()
        if not items and self.file_list.currentItem() is not None:
            items = [self.file_list.currentItem()]
        return [Path(item.data(Qt.ItemDataRole.UserRole)) for item in items]

    def _open_selected_files(self) -> None:
        paths = self._selected_file_paths()
        if not paths:
            self._warn("Select at least one FITS file first.")
            return

        for path in paths:
            loaded = iter_displayable_slices(path)
            self.loaded_files[path] = loaded
            self._update_file_item_status(path, loaded.status)
            if loaded.status == "ok":
                self._log(f"Opened {path.name}: {len(loaded.slices)} displayable slice(s).")
            else:
                self._log(f"{path.name}: {loaded.status}: {loaded.message}")
        self._refresh_slice_list()

    def _update_file_item_status(self, path: Path, status: str) -> None:
        for row in range(self.file_list.count()):
            item = self.file_list.item(row)
            if Path(item.data(Qt.ItemDataRole.UserRole)) == path:
                prefix = {"ok": "[ok]", "invalid": "[invalid]", "no-image": "[no image]"}.get(status, "[?]")
                item.setText(f"{prefix} {path.name}")
                break

    def _refresh_slice_list(self) -> None:
        self.slice_list.clear()
        for loaded in self.loaded_files.values():
            if loaded.status != "ok":
                continue
            for image_slice in loaded.slices:
                item = QListWidgetItem(f"{image_slice.file_path.name}: {image_slice.label}")
                item.setData(Qt.ItemDataRole.UserRole, image_slice)
                self.slice_list.addItem(item)

    def _discard_selected_file(self) -> None:
        current = self.file_list.currentItem()
        if current is None:
            return
        path = Path(current.data(Qt.ItemDataRole.UserRole))
        self.loaded_files.pop(path, None)
        self.file_list.takeItem(self.file_list.row(current))
        self._refresh_slice_list()
        self._log(f"Discarded {path.name} from this EDA session.")

    def _select_slice(self) -> None:
        item = self.slice_list.currentItem()
        if item is None:
            return
        image_slice = item.data(Qt.ItemDataRole.UserRole)
        self.current_slice = image_slice
        self.current_image_adu = image_slice.data.copy()
        self.current_image_electrons = None
        self.background_result = None
        self.psf_result = None
        self._set_editing_enabled(False)
        self.bg_label.setText("No background estimate yet.")
        self.psf_label.setText("No PSF estimate yet.")

        ny, nx = self.current_image_adu.shape
        self.center_x_spin.setRange(0, nx - 1)
        self.center_y_spin.setRange(0, ny - 1)
        self.center_x_spin.setValue(nx // 2)
        self.center_y_spin.setValue(ny // 2)
        self.crop_nx_spin.setMaximum(nx)
        self.crop_ny_spin.setMaximum(ny)
        self.canvas.show_image(self.current_image_adu, f"{image_slice.file_path.name}: {image_slice.label}")
        self._log(f"Selected {image_slice.file_path.name}, {image_slice.label}, shape={self.current_image_adu.shape}.")

    def _start_editing(self) -> None:
        if self.current_slice is None or self.current_image_adu is None:
            self._warn("Select a displayable FITS slice first.")
            return
        self._set_editing_enabled(True)
        self._preview_crop()
        self._log("Editing enabled. Click the image to set crop center.")

    def _compute_background(self) -> None:
        if self.current_slice is None or self.current_image_adu is None:
            self._warn("Select an image slice first.")
            return
        try:
            result = estimate_background_and_gain(self.current_image_adu, self.current_slice.header)
        except Exception as exc:
            self._warn(f"Background/gain estimation failed: {exc}")
            return

        self.background_result = result
        self.current_image_electrons = result.image_electrons
        self.canvas.show_image(self.current_image_electrons, f"{self.current_slice.file_path.name}: electrons")
        self._preview_crop()
        self.bg_label.setText(
            "Background ADU mean={:.6g}, sigma={:.6g}\n"
            "Gain={:.6g} ({})\n"
            "Electron mean={:.6g}, sigma={:.6g}\n"
            "Masked bright/object fraction={:.3f}".format(
                result.background_mean_adu,
                result.background_sigma_adu,
                result.gain,
                result.gain_source,
                result.background_mean_electrons,
                result.background_sigma_electrons,
                result.masked_fraction,
            )
        )
        self._log("Computed background/gain and updated display to electron counts.")

    def _compute_psf(self) -> None:
        image = self.current_image_electrons if self.current_image_electrons is not None else self.current_image_adu
        if image is None:
            self._warn("Select an image slice first.")
            return
        try:
            result = estimate_psf(image, self.pixel_scale_spin.value())
        except Exception as exc:
            self._warn(f"PSF estimation failed: {exc}")
            return

        self.psf_result = result
        self.psf_label.setText(
            "Sources detected={}\n"
            "Gaussian fits used={}\n"
            "PSF sigma={:.6g} pix\n"
            "Pipeline PSF FWHM={:.6g} arcsec".format(
                result.n_sources,
                result.n_fitted,
                result.sigma_pix,
                result.fwhm_arcsec,
            )
        )
        self._log("Estimated PSF from detected stars.")

    def _set_center_from_click(self, x: int, y: int) -> None:
        if self.current_image_adu is None:
            return
        ny, nx = self.current_image_adu.shape
        self.center_x_spin.setValue(int(np.clip(x, 0, nx - 1)))
        self.center_y_spin.setValue(int(np.clip(y, 0, ny - 1)))
        self._preview_crop()

    def _preview_crop(self) -> None:
        if self.current_image_adu is None:
            return
        self.canvas.set_crop_rect(
            self.center_x_spin.value(),
            self.center_y_spin.value(),
            self.crop_nx_spin.value(),
            self.crop_ny_spin.value(),
        )

    def _crop_bounds(self) -> tuple[int, int, int, int]:
        if self.current_image_adu is None:
            raise ValueError("No current image.")
        ny_image, nx_image = self.current_image_adu.shape
        nx = self.crop_nx_spin.value()
        ny = self.crop_ny_spin.value()
        cx = self.center_x_spin.value()
        cy = self.center_y_spin.value()
        x0 = int(round(cx - nx / 2))
        y0 = int(round(cy - ny / 2))
        x0 = max(0, min(x0, nx_image - nx))
        y0 = max(0, min(y0, ny_image - ny))
        x1 = x0 + nx
        y1 = y0 + ny
        if x0 < 0 or y0 < 0 or x1 > nx_image or y1 > ny_image:
            raise ValueError("Crop is larger than the image.")
        return x0, x1, y0, y1

    def _save_crop(self) -> None:
        if self.current_slice is None:
            self._warn("Select and edit an image slice before saving.")
            return
        if self.background_result is None or self.current_image_electrons is None:
            self._warn("Compute background/gain before saving so the crop is in electron counts.")
            return
        if self.psf_result is None:
            self._warn("Estimate PSF before saving so the manifest has PSF parameters.")
            return

        nx = self.crop_nx_spin.value()
        ny = self.crop_ny_spin.value()
        requested_shape = (ny, nx)
        if self.locked_crop_shape is not None and requested_shape != self.locked_crop_shape and not self.allow_shape_reset.isChecked():
            self._warn(
                f"This output project is locked to shape {self.locked_crop_shape}. "
                "Choose a new output folder or enable Allow new project shape."
            )
            return

        output_dir = Path(self.output_edit.text()).expanduser()
        if not output_dir.is_absolute():
            output_dir = SCRIPT_DIR / output_dir
        self.output_dir = output_dir
        images_dir = self.output_dir / "synthetic" / "images"
        images_dir.mkdir(parents=True, exist_ok=True)

        x0, x1, y0, y1 = self._crop_bounds()
        sample_id = self._unique_sample_id(safe_sample_id(self.current_slice.file_path, self.current_slice.label))
        npy_path = images_dir / f"{sample_id}.npy"
        try:
            if self._manifest_path().exists():
                self._validate_manifest_shape(self._manifest_path(), ny, nx)
        except Exception as exc:
            self._warn(f"Cannot append to manifest: {exc}")
            return

        crop = self.current_image_electrons[y0:y1, x0:x1].astype(np.float32)
        np.save(npy_path, crop)

        row = {
            "sample_id": sample_id,
            "image_path": str(npy_path),
            "truth_json_path": "",
            "env_nx": nx,
            "env_ny": ny,
            "env_pixel_scale": self.pixel_scale_spin.value(),
            "env_noise_sigma": self.background_result.background_sigma_electrons,
            "env_psf_sigma": self.psf_result.fwhm_arcsec,
        }
        self._append_manifest_row(row)
        self.locked_crop_shape = requested_shape
        self.allow_shape_reset.setChecked(False)
        self._update_fit_command()
        self._log(f"Saved crop {sample_id}: {crop.shape} -> {npy_path}")

    def _unique_sample_id(self, base: str) -> str:
        manifest_path = self._manifest_path()
        existing: set[str] = set()
        if manifest_path.exists():
            with manifest_path.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle)
                existing = {row["sample_id"] for row in reader if row.get("sample_id")}
        candidate = base
        idx = 1
        while candidate in existing or (self.output_dir / "synthetic" / "images" / f"{candidate}.npy").exists():
            candidate = f"{base}_{idx:03d}"
            idx += 1
        return candidate

    def _append_manifest_row(self, row: dict[str, Any]) -> None:
        manifest_path = self._manifest_path()
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        write_header = not manifest_path.exists()
        if not write_header:
            self._validate_manifest_shape(manifest_path, int(row["env_ny"]), int(row["env_nx"]))
        with manifest_path.open("a", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS)
            if write_header:
                writer.writeheader()
            writer.writerow(row)

    def _validate_manifest_shape(self, manifest_path: Path, ny: int, nx: int) -> None:
        with manifest_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                if not row:
                    continue
                prev_shape = (int(float(row["env_ny"])), int(float(row["env_nx"])))
                if prev_shape != (ny, nx) and not self.allow_shape_reset.isChecked():
                    raise ValueError(f"Existing manifest shape is {prev_shape}, requested {(ny, nx)}.")
                break

    def _manifest_path(self) -> Path:
        output_dir = Path(self.output_edit.text()).expanduser()
        if not output_dir.is_absolute():
            output_dir = SCRIPT_DIR / output_dir
        return output_dir / "synthetic_manifest.csv"

    def _read_locked_crop_shape(self) -> Optional[tuple[int, int]]:
        manifest_path = DEFAULT_OUTPUT_DIR / "synthetic_manifest.csv"
        if not manifest_path.exists():
            return None
        try:
            with manifest_path.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle)
                row = next(reader, None)
                if row is None:
                    return None
                return int(float(row["env_ny"])), int(float(row["env_nx"]))
        except Exception:
            return None

    def _choose_output_dir(self) -> None:
        directory = QFileDialog.getExistingDirectory(self, "Choose output folder", str(self.output_dir))
        if not directory:
            return
        self.output_dir = Path(directory)
        self.output_edit.setText(str(self.output_dir))
        self.locked_crop_shape = self._read_shape_for_output(self.output_dir)
        self._update_fit_command()

    def _read_shape_for_output(self, output_dir: Path) -> Optional[tuple[int, int]]:
        manifest_path = output_dir / "synthetic_manifest.csv"
        if not manifest_path.exists():
            return None
        try:
            with manifest_path.open("r", encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle), None)
            if row is None:
                return None
            return int(float(row["env_ny"])), int(float(row["env_nx"]))
        except Exception:
            return None

    def _update_fit_command(self) -> None:
        output_dir = Path(self.output_edit.text()).expanduser()
        output_arg_raw = output_dir.name if output_dir.is_absolute() and output_dir.parent == SCRIPT_DIR else str(output_dir)
        output_arg = shlex.quote(output_arg_raw)
        command = (
            "python fit_quality_pipeline.py --mode fit "
            f"--output-dir {output_arg} "
            f"--nx {self.crop_nx_spin.value()} "
            f"--ny {self.crop_ny_spin.value()} "
            f"--pixel-scale {self.pixel_scale_spin.value():.5g}"
        )
        self.fit_command_label.setText(command + "\n\nUse --mode fit, not --mode all, for exported real data.")

    def _log(self, message: str) -> None:
        self.status_box.append(message)

    def _warn(self, message: str) -> None:
        self._log(message)
        QMessageBox.warning(self, "Photometry EDA", message)


def main() -> None:
    app = QApplication(sys.argv)
    window = PhotometryEdaWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
