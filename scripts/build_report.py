"""Build the submission report PDF, with six vector figure pages and their captions."""
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
'Fig3':'Separate discovery, IFN-programme and protein measurements',
'Fig4':'Candidate regulators and expression-associated genetic evidence',
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
 if need>iy-58:raise RuntimeError('Caption does not fit on the page: '+name)
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
with (out/'Nature_figures_with_legends.pdf').open('wb') as f:writer.write(f)
# Build the narrative PDF directly, so that table text stays inside its cells.
from reportlab.platypus import SimpleDocTemplate, Spacer, Table, TableStyle, PageBreak, KeepTogether
from reportlab.lib import colors
import re
story=(root/'docs/story.md').read_text().split('## Figure legends')[0]
validation=(root/'docs/REVISION_VALIDATION.md').read_text()
body=story+'\n\n## Final classifier and validation records\n\nSee classifier_baselines.tsv and classifier_protocol.json for all estimates and the resampling protocol.\n\n'+validation
styles={
 'body':ParagraphStyle('body',fontName='ReportArial',fontSize=9.7,leading=13.2,spaceAfter=6),
 'h1':ParagraphStyle('h1',fontName='ReportArialBold',fontSize=19,leading=23,textColor=colors.HexColor('#202020'),spaceBefore=10,spaceAfter=10,keepWithNext=True),
 'h2':ParagraphStyle('h2',fontName='ReportArialBold',fontSize=13.5,leading=17,textColor=colors.HexColor('#202020'),spaceBefore=10,spaceAfter=7,keepWithNext=True),
 'h3':ParagraphStyle('h3',fontName='ReportArialBold',fontSize=11,leading=14,spaceBefore=7,spaceAfter=5,keepWithNext=True),
 'cell':ParagraphStyle('cell',fontName='ReportArial',fontSize=8.3,leading=10.8),
 'head':ParagraphStyle('head',fontName='ReportArialBold',fontSize=8.3,leading=10.8,textColor=colors.HexColor('#202020')),
}
def inline(text):
 text=html.escape(text,quote=False)
 text=re.sub(r'\[([^\]]+)\]\((https?://[^)]+)\)',r'<link href="\2" color="#202020">\1</link>',text)
 text=re.sub(r'\*\*(.+?)\*\*',r'<b>\1</b>',text)
 text=re.sub(r'`([^`]+)`',r'\1',text)
 return text
flow=[];lines=body.splitlines();i=0
while i<len(lines):
 line=lines[i].strip()
 if not line:i+=1;continue
 if line.startswith('|'):
  rows=[]
  while i<len(lines) and lines[i].strip().startswith('|'):
   cells=[c.strip() for c in lines[i].strip().strip('|').split('|')]
   if not all(re.fullmatch(r':?-+:?',c) for c in cells):rows.append(cells)
   i+=1
  n=len(rows[0]);widths=[108,225,186] if n==3 else [115,280,65,59]
  data=[[Paragraph(inline(c),styles['head' if j==0 else 'cell']) for c in row] for j,row in enumerate(rows)]
  table=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT')
  table.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#EEEEEE')),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEABOVE',(0,0),(-1,0),.65,colors.HexColor('#444444')),('LINEBELOW',(0,0),(-1,0),.5,colors.HexColor('#888888')),('LINEBELOW',(0,-1),(-1,-1),.65,colors.HexColor('#444444')),('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.white])]))
  flow.extend([table,Spacer(1,8)]);continue
 if line.startswith('#'):
  level=len(line)-len(line.lstrip('#'));text=line[level:].strip()
  if text=='Final classifier and validation records':flow.append(PageBreak())
  flow.append(Paragraph(inline(text),styles['h'+str(min(level,3))]));i+=1;continue
 if line.startswith('- '):
  flow.append(Paragraph('• '+inline(line[2:]),styles['body']));i+=1;continue
 paragraph=[line];i+=1
 while i<len(lines) and lines[i].strip() and not lines[i].lstrip().startswith(('#','|','- ')):
  paragraph.append(lines[i].strip());i+=1
 flow.append(Paragraph(inline(' '.join(paragraph)),styles['body']))
def footer(c,doc):
 c.setFont('ReportArial',7.5);c.setFillColor(colors.HexColor('#606060'))
 c.drawString(38,25,'GeNeDis 2026 · Xia Wu · University of Iowa')
 c.drawRightString(W-38,25,str(doc.page))
narrative=tmp/'narrative.pdf'
SimpleDocTemplate(str(narrative),pagesize=A4,leftMargin=38,rightMargin=38,topMargin=36,bottomMargin=40,title='GeNeDis 2026: RA multi-omics submission',author='Xia Wu').build(flow,onFirstPage=footer,onLaterPages=footer)
combined=PdfWriter()
for path in [narrative,out/'Nature_figures_with_legends.pdf']:
 for page in PdfReader(path).pages:combined.add_page(page)
with (out/'GeNeDis2026_submission_report.pdf').open('wb') as f:combined.write(f)
# The report is submitted as a PDF, so remove any out-of-date editable copy.
(out/'GeNeDis2026_submission_report.docx').unlink(missing_ok=True)
print('Saved submission report:',len(combined.pages),'pages, including six vector figure pages with captions')
