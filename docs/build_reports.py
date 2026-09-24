"""Сборка всех Markdown-отчётов в PDF с формулами mathtext."""
from __future__ import annotations

import argparse
import base64
import io
import re
from pathlib import Path

import markdown
from matplotlib.mathtext import math_to_image
from weasyprint import HTML

CSS = '''
@page { size: A4; margin: 19mm 18mm 18mm 18mm; @bottom-center { content: counter(page); font-size: 9pt; color: #64748b; } }
body { font-family: "DejaVu Sans", sans-serif; color: #162033; font-size: 10pt; line-height: 1.48; }
h1 { font-size: 20pt; color: #173455; margin-bottom: 9mm; border-bottom: 1.5pt solid #547999; padding-bottom: 3mm; }
h2 { font-size: 14pt; color: #20476a; margin-top: 8mm; break-after: avoid; }
h3 { font-size: 11.5pt; color: #315a79; break-after: avoid; }
p { orphans: 3; widows: 3; }
img { max-width: 100%; height: auto; }
.equation { text-align: center; margin: 4mm 0; }
.equation img { max-height: 23mm; }
.math-inline { vertical-align: middle; max-height: 13pt; }
table { width: 100%; border-collapse: collapse; font-size: 8.3pt; margin: 4mm 0; }
th,td { border: .5pt solid #b9c8d5; padding: 2mm; vertical-align: top; overflow-wrap: anywhere; }
th { background: #e9f1f7; }
pre,code { font-family: "DejaVu Sans Mono", monospace; font-size: 8pt; }
pre { white-space: pre-wrap; background: #f2f5f7; padding: 3mm; }
blockquote { border-left: 2pt solid #7ba4bf; padding-left: 4mm; color: #405365; }
a { color: #245b88; text-decoration: underline; }
figure, img { break-inside: avoid; }
tr { break-inside: avoid; }
table td:first-child, table th:first-child { min-width: 72pt; }
'''


def formula_image(formula: str, display: bool) -> str:
    stream=io.BytesIO()
    math_to_image(f'${formula.strip()}$',stream,dpi=180,format='svg',color='#162033')
    uri='data:image/svg+xml;base64,'+base64.b64encode(stream.getvalue()).decode('ascii')
    return (f'<div class="equation"><img src="{uri}" alt="{formula}"/></div>' if display
            else f'<img class="math-inline" src="{uri}" alt="{formula}"/>')


def convert(md: str) -> str:
    expressions=[]
    def block(match):
        expressions.append(formula_image(match.group(1),True))
        return f'@@FORMULA{len(expressions)-1}@@'
    def inline(match):
        expressions.append(formula_image(match.group(1),False))
        return f'@@FORMULA{len(expressions)-1}@@'
    md=re.sub(r'\$\$\s*(.*?)\s*\$\$',block,md,flags=re.S)
    md=re.sub(r'(?<!\$)\$([^$\n]+)\$(?!\$)',inline,md)
    html=markdown.markdown(md,extensions=['tables','fenced_code','toc'])
    for i,expr in enumerate(expressions):html=html.replace(f'@@FORMULA{i}@@',expr)
    return '<!doctype html><html lang="ru"><meta charset="utf-8"><style>'+CSS+'</style><body>'+html+'</body></html>'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('paths',nargs='*',type=Path)
    args=p.parse_args()
    paths=args.paths or sorted((Path(__file__).parent/'reports').glob('*.md'))
    for path in paths:
        html=convert(path.read_text(encoding='utf-8'))
        output=path.with_suffix('.pdf')
        HTML(string=html,base_url=str(path.parent.resolve())).write_pdf(output)
        print(f'[pdf] {output} ({output.stat().st_size} bytes)')


if __name__=='__main__':main()
