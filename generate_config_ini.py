import glob
import re
import sys
import datetime

try:
    indir = sys.argv[1]
except:
    indir = '.'
folders = glob.glob(indir + '/*/')
folders = [x for x in folders if not re.search('backup', x)]
bedgraphs = [x for x in folders if (
    re.search('bedgraph', x) or re.search('wig', x))]
#norm_bedgraph = [x for x in bedgraphs if (
#    re.search('norm', x) and not re.search('unnorm', x))]
unnorm_bedgraph = [x for x in bedgraphs if (
    re.search('unnorm', x))]
bed = [x for x in folders if (
    re.search('bed', x) and (x not in bedgraphs) and not re.search('uncollapse', x))]
for fn in [bed,  unnorm_bedgraph]:
    if len(fn) > 1:
        print fn
        sys.exit()
bed = bed[0]
#norm_bedgraph = norm_bedgraph[0]
unnorm_bedgraph = unnorm_bedgraph[0]
print folders
print bedgraphs
#print norm_bedgraph
print unnorm_bedgraph
print bed
controlbeds = [x for x in glob.glob(bed + '/*.bed') if re.search('control', x)]
controlbeds = sorted(controlbeds, key=lambda x: len(x))
rnaseqbeds = [x for x in glob.glob(bed + '/*.bed') if re.search('seq', x)]
rnaseqbeds = sorted(rnaseqbeds, key=lambda x: len(x))
exp_beds = [x for x in glob.glob(bed + '/*.bed') if re.search('exp', x)]
clip_reps = [x for x in glob.glob(unnorm_bedgraph + '/*.wig') if \
    re.search('exp', x)]
clip_reps = list(set([re.sub('_[\+-]\.wig', '.wig', x) for x in \
                 clip_reps]))
#for rep in clip_reps:
#    rep = re.sub('_[\+-].wig', '.wig', rep)
li = """
[library]
top: {top}
lib: {lib}
gtf: {gtf}
gtf_filename: {gtf}
bedgraphs_folder: {unnorm}
read_beds: {collapsed}
bed_dir: {collapsed}
control_bed1: {controlbed}
exp_bed1: {e1bed}
exp_bed2: {e2bed}
exp_bed3: {e3bed}
experiment_name: {exp_name}
min_rep_number: 2
clip_replicate1: {cr1}
clip_replicate2: {cr2}
clip_replicate3: {cr3}
neg_ip_filename: {controlbed}
rna_seq_filename: {seq}
positive_control_genes: gld-1,htp-1,htp-2,mpk-1,him-3,fbf-1,lip-1,syp-2,fbf-2,fog-1,fem-3,syp-3,gld-3,fog-3,egl-4
motif: tgt\w\w\wat
""".format(
top=indir,
lib='/groups/Kimble/Common/fog_iCLIP/calls/lib/',
gtf="%(lib)s/gtf_with_names_column.txt",
unnorm=unnorm_bedgraph,
collapsed=bed,
controlbed=controlbeds[0],
e1bed=exp_beds[0],
e2bed=exp_beds[1],
e3bed=exp_beds[2],
exp_name=datetime.datetime.now().strftime("%I_%M_%p"),
cr1=clip_reps[0],
cr2=clip_reps[1],
cr3=clip_reps[2],
seq=rnaseqbeds[0]
)
print li
with open('auto.ini', 'w') as f: f.write(li)
