"""Построить условные расстояния по линейному закону Хаббла из проверенной таблицы z."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

C_KM_S = 299792.458
ROOT = Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--redshifts', type=Path, default=ROOT/'quasar_redshifts.csv')
    p.add_argument('--output', type=Path, default=ROOT/'quasar_distances_hubble.csv')
    p.add_argument('--h0', type=float, default=70.0, help='км/с/Мпк; для отчёта используется 70')
    args = p.parse_args()
    if args.h0 <= 0: p.error('H0 должен быть положительным')
    with args.redshifts.open(encoding='utf-8',newline='') as f:
        sources = list(csv.DictReader(f))
    with args.output.open('w',encoding='utf-8',newline='') as f:
        writer=csv.writer(f,lineterminator='\n')
        writer.writerow(['object','z','h0_km_s_mpc','hubble_distance_mpc',
                         'hubble_distance_gpc','source_url','match_note'])
        for row in sources:
            z=float(row['z'])
            d=C_KM_S*z/args.h0
            writer.writerow([row['object'],row['z'],f'{args.h0:g}',f'{d:.2f}',
                             f'{d/1000:.5f}',row['source_url'],row['match_note']])
    print(f'{len(sources)} источника: {args.output}')


if __name__=='__main__':main()
