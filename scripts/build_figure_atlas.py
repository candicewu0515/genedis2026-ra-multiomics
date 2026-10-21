"""Build the separate six-page figure atlas from accepted vector figures."""
from pathlib import Path
import sys,json,subprocess,tempfile,html,os,shutil
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from pypdf import PdfReader,PdfWriter,Transformation
root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parents[1]
out=root/'presentation';out.mkdir(exist_ok=True)
fig=root/'results/figures/main'
font=Path('/System/Library/Fonts/Supplemental/Arial.ttf')
bold=Path('/System/Library/Fonts/Supplemental/Arial Bold.ttf')
if not font.exists():
 from matplotlib.font_manager import findfont
 font=Path(findfont('DejaVu Sans'));bold=Path(findfont('DejaVu Sans:bold'))
pdfmetrics.registerFont(TTFont('ReportArial',str(font)))
pdfmetrics.registerFont(TTFont('ReportArialBold',str(bold)))
pdfmetrics.registerFontFamily('ReportArial',normal='ReportArial',bold='ReportArialBold',italic='ReportArial',boldItalic='ReportArialBold')
tmp=Path(tempfile.mkdtemp(prefix='genedis_report_'))
guide=json.loads((fig/'panel_guide.json').read_text())
titles={'Fig1':'Study design, published atlas states and integration baselines',
'Fig2':'Genetic association and spatial neighbourhoods',
'Fig3':'Distinct discovery, IFN-programme and protein measurements',
'Fig4':'Candidate regulatory and expression-associated genetic evidence',
'ExtData_Fig1':'Baselines, robustness and boundary findings',
'ExtData_Fig2':'Exploratory synovial-fluid and NSCLC extensions'}
writer=PdfWriter();W,H=A4;left=38;top=H-42;width=W-2*left
for name,panels in guide.items():
 original=PdfReader(fig/(name+'.pdf')).pages[0]
 scale=width/float(original.mediabox.width);ih=float(original.mediabox.height)*scale
 iy=top-ih
 caption=[f'<b>{html.escape(k)}</b>  {html.escape(v)}' for k,v in panels.items()]
 fontsize=8.3
 while True:
  style=ParagraphStyle('caption',fontName='ReportArial',fontSize=fontsize,leading=fontsize+1.8,spaceAfter=3.1)
  parts=[Paragraph(x,style) for x in caption]
  sizes=[p.wrap(width,H) for p in parts]
  need=sum(h+style.spaceAfter for _,h in sizes)
  if need <= iy-58 or fontsize<=7.2:break
  fontsize-=.2
 if need>iy-58:raise RuntimeError('Caption exceeds available page: '+name)
 base=tmp/(name+'_caption.pdf');c=canvas.Canvas(str(base),pagesize=A4)
 c.setFont('ReportArialBold',10.8);c.setFillColorRGB(.12,.12,.12)
 c.drawString(left,H-27,name.replace('_',' ')+' | '+titles[name])
 y=iy-9
 for p,(_,h) in zip(parts,sizes):
  p.drawOn(c,left,y-h);y-=h+style.spaceAfter
 c.setFont('ReportArial',7.5);c.setFillColorRGB(.35,.35,.35)
 c.drawString(left,27,'GeNeDis 2026 · Xia Wu · University of Iowa')
 c.save()
 page=PdfReader(base).pages[0];page.merge_transformed_page(original,Transformation().scale(scale).translate(left,iy),over=True)
 writer.add_page(page)
with (out/'GeNeDis2026_figure_atlas.pdf').open('wb') as f:writer.write(f)

print('Saved six-page figure atlas')
