"""Build the competition technical report and a separate Methods/QC appendix.

Narrative sources are separate from the historical analysis/revision records.
Scientific figures are supplied in the separate, accepted six-page atlas.
"""
from pathlib import Path
import sys, html, re, json, subprocess, hashlib
import pandas as pd
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, Flowable
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from pypdf import PdfReader
ROOT=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parents[1]
OUT=ROOT/'presentation';OUT.mkdir(exist_ok=True)
font=Path('/System/Library/Fonts/Supplemental/Arial.ttf');bold=Path('/System/Library/Fonts/Supplemental/Arial Bold.ttf')
if not font.exists():
 from matplotlib.font_manager import findfont
 font=Path(findfont('DejaVu Sans'));bold=Path(findfont('DejaVu Sans:bold'))
pdfmetrics.registerFont(TTFont('ReportArial',str(font)));pdfmetrics.registerFont(TTFont('ReportArialBold',str(bold)))
pdfmetrics.registerFontFamily('ReportArial',normal='ReportArial',bold='ReportArialBold',italic='ReportArial',boldItalic='ReportArialBold')
W,H=A4;M=44;WIDTH=W-2*M
INK=colors.HexColor('#202020');GREY=colors.HexColor('#606060')
ST={
 'body':ParagraphStyle('body',fontName='ReportArial',fontSize=10.25,leading=13.9,spaceAfter=6,allowWidows=0,allowOrphans=0),
 'meta':ParagraphStyle('meta',fontName='ReportArial',fontSize=9.5,leading=12.8,spaceAfter=6),
 'h1':ParagraphStyle('h1',fontName='ReportArialBold',fontSize=20,leading=24,spaceAfter=10,keepWithNext=True,textColor=INK),
 'h2':ParagraphStyle('h2',fontName='ReportArialBold',fontSize=14,leading=18,spaceAfter=8,spaceBefore=3,keepWithNext=True,textColor=INK),
 'h3':ParagraphStyle('h3',fontName='ReportArialBold',fontSize=11.2,leading=14.3,spaceAfter=5,spaceBefore=6,keepWithNext=True,textColor=INK),
 'cell':ParagraphStyle('cell',fontName='ReportArial',fontSize=9.0,leading=11.8),
 'head':ParagraphStyle('head',fontName='ReportArialBold',fontSize=9.0,leading=11.8,textColor=INK),
 'ref':ParagraphStyle('ref',fontName='ReportArial',fontSize=8.6,leading=11.3,spaceAfter=4),
 'bullet':ParagraphStyle('bullet',fontName='ReportArial',fontSize=10.25,leading=13.9,spaceAfter=5,leftIndent=9,firstLineIndent=-8),
}
def inline(text):
 text=html.escape(text,quote=False)
 text=re.sub(r'\[([^\]]+)\]\((https?://[^)]+)\)',r'<link href="\2" color="#202020">\1</link>',text)
 text=re.sub(r'\*\*(.+?)\*\*',r'<b>\1</b>',text)
 return re.sub(r'`([^`]+)`',r'\1',text)
class Workflow(Flowable):
 def __init__(self):super().__init__();self.width=WIDTH;self.height=80
 def draw(self):
  c=self.canv;gap=13;bw=(WIDTH-2*gap)/3;bh=61
  items=[('Open resources','RNA/ATAC · bulk · spatial','Protein · GWAS/QTL'),('Cohort-specific analysis','QC and measured features','Integration and model baselines'),('Evidence and delivery','Tables and source-linked results','Figure atlas and follow-up hypotheses')]
  for i,(title,a,b) in enumerate(items):
   x=i*(bw+gap);c.setFillColor(colors.HexColor('#F0F0F0'));c.setStrokeColor(colors.HexColor('#AAAAAA'));c.setLineWidth(.45);c.roundRect(x,12,bw,bh,4,fill=1,stroke=1)
   c.setFillColor(INK);c.setFont('ReportArialBold',9.0);c.drawCentredString(x+bw/2,55,title)
   c.setFont('ReportArial',8.1);c.drawCentredString(x+bw/2,39,a);c.drawCentredString(x+bw/2,26,b)
   if i<2:
    c.setStrokeColor(GREY);c.line(x+bw+2,42,x+bw+gap-2,42);c.line(x+bw+gap-5,45,x+bw+gap-2,42);c.line(x+bw+gap-5,39,x+bw+gap-2,42)

def classifier_table():
 d=pd.read_csv(ROOT/'results/tables/classifier_baselines.tsv',sep='\t');tasks=['T1 RA vs OA','T2 ACPA+ vs ACPA- (RA)'];models=['F0 covariates','F1 eight programme scores','F2 top-30 prioritised genes','F3 genome-wide elastic net']
 rows=['| Model | RA vs OA | ACPA+ vs ACPA- |','|---|---|---|']
 for model in models:
  values=[]
  for task in tasks:
   x=d[(d.task==task)&(d.model==model)].iloc[0];values.append(f'{x.auc_mean:.3f} ± {x.auc_sd:.3f}')
  rows.append('| '+model+' | '+' | '.join(values)+' |')
 return '\n'.join(rows)

def table_flow(rows):
 head=rows[0][0];n=len(rows[0])
 if head=='Model':widths=[.46,.27,.27]
 elif head=='Resource':widths=[.20,.20,.26,.34]
 elif head=='Measurement':widths=[.20,.39,.41]
 elif head=='Component':widths=[.19,.36,.45]
 elif head=='Score':widths=[.14,.27,.26,.33]
 elif head=='Output':widths=[.30,.23,.47]
 else:widths=[1/n]*n
 data=[[Paragraph(inline(x),ST['head' if j==0 else 'cell']) for x in row] for j,row in enumerate(rows)]
 t=Table(data,colWidths=[WIDTH*x for x in widths],repeatRows=1,hAlign='LEFT')
 t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#EEEEEE')),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEABOVE',(0,0),(-1,0),.65,colors.HexColor('#444444')),('LINEBELOW',(0,0),(-1,0),.45,colors.HexColor('#888888')),('LINEBELOW',(0,-1),(-1,-1),.65,colors.HexColor('#444444')),('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]))
 return t

def parse(text):
 lines=text.splitlines();flow=[];i=0;meta=True
 while i<len(lines):
  line=lines[i].strip()
  if not line:i+=1;continue
  if line=='<!-- PAGEBREAK -->':flow.append(PageBreak());i+=1;meta=False;continue
  if line=='{{WORKFLOW}}':flow.append(Workflow());i+=1;continue
  if line.startswith('|'):
   rows=[]
   while i<len(lines) and lines[i].strip().startswith('|'):
    cells=[c.strip() for c in lines[i].strip().strip('|').split('|')]
    if not all(re.fullmatch(r':?-+:?',c) for c in cells):rows.append(cells)
    i+=1
   flow.extend([table_flow(rows),Spacer(1,8)]);continue
  if line.startswith('#'):
   level=len(line)-len(line.lstrip('#'));flow.append(Paragraph(inline(line[level:].strip()),ST['h'+str(min(level,3))]));meta=(level==1);i+=1;continue
  if line.startswith('- '):flow.append(Paragraph('• '+inline(line[2:]),ST['bullet']));i+=1;continue
  p=[line];i+=1
  while i<len(lines) and lines[i].strip() and not lines[i].lstrip().startswith(('#','|','- ','<!--','{{')):p.append(lines[i].strip());i+=1
  text=' '.join(p);sty='ref' if re.match(r'^\[\d+\]',text) else ('meta' if meta else 'body')
  flow.append(Paragraph(inline(text),ST[sty]))
 return flow

def build(source,name,label):
 text=(ROOT/'docs'/source).read_text().replace('{{CLASSIFIER_TABLE}}',classifier_table())
 def page(c,doc):
  c.setFont('ReportArial',7.5);c.setFillColor(GREY);c.drawString(M,24,'GeNeDis 2026 · Xia Wu · University of Iowa');c.drawRightString(W-M,24,label+' | '+str(doc.page))
  if doc.page>1:c.drawString(M,H-25,'RA multi-omics workflow')
 target=OUT/name
 SimpleDocTemplate(str(target),pagesize=A4,leftMargin=M,rightMargin=M,topMargin=41,bottomMargin=39,title='GeNeDis 2026: '+label,author='Xia Wu').build(parse(text),onFirstPage=page,onLaterPages=page)
 return target

atlas=OUT/'GeNeDis2026_figure_atlas.pdf'
if not atlas.exists():subprocess.run([sys.executable,str(ROOT/'scripts/build_figure_atlas.py'),str(ROOT)],check=True)
assert len(PdfReader(atlas).pages)==6
report=build('competition_report.md','GeNeDis2026_submission_report.pdf','Technical report')
methods=build('methods_and_QC.md','GeNeDis2026_methods_QC.pdf','Methods and QC')
check={}
for name,p in [('report',report),('methods_QC',methods)]:
 reader=PdfReader(p);text='\n'.join(x.extract_text() for x in reader.pages)
 assert 'xwu76@uiowa.edu' in text and '{{' not in text and '<!--' not in text
 assert 'Revision traceability' not in text and 'early step-36 branch' not in text
 for page in reader.pages:assert abs(float(page.mediabox.width)-W)<1 and abs(float(page.mediabox.height)-H)<1
 check[name]={'pages':len(reader.pages),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
print(json.dumps(check,indent=2))
(ROOT/'results/report_document_checks.json').write_text(json.dumps({'status':'PASS','scope':'document generation, source placeholders and page format; interpretation reviewed separately','documents':check},indent=2))
