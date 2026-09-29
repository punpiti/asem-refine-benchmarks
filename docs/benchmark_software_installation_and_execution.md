# การติดตั้งและรันซอฟต์แวร์สำหรับ benchmark

เอกสารนี้สรุปชุด benchmark ของ manuscript ASEM/ASEM-Hybrid ว่าโปรแกรมมาจากไหน
ติดตั้งไว้ที่ใด wrapper เรียกอย่างไร ตั้ง parameter อะไร และเก็บผลลัพธ์ที่ไหน
ข้อมูลตรวจจาก source code, environment history และ package metadata ที่อยู่ในเครื่อง
ณ วันที่ 29 กันยายน 2026

## 1. ภาพรวม

repository ของ benchmark อยู่ที่:

```text
/home/punpiti/OneDrive/genome/asem-refine-benchmarks
```

| วิธี | รุ่นที่ใช้จริง | environment/ตำแหน่ง |
|---|---:|---|
| ASEM non-recursive | project source | `scripts/baseline_ieee_access.py` |
| ASEM recursive | project source | `scripts/baseline_ojemb.py` |
| SRSC | project source | `scripts/baseline_ecticon.py` |
| ASEM-Hybrid | project source | `scripts/asem_hybrid.py` และ boundary variant |
| NOVOPlasty | 4.3.5 | `~/.local/share/mamba/envs/genome/` |
| GetOrganelle | 1.7.7.0 | `~/.local/share/mamba/envs/getorganelle/` |
| MIA | 1.0 | `~/.local/share/mamba/envs/mia-assembler/` |
| Pilon | 1.24 | `~/.local/share/mamba/envs/pilon-benchmark/` |
| MITObim | 1.9.1 | `~/.local/share/mamba/envs/mitobim/`; ติดตั้งได้แต่ไม่รวมในผล |

วิธีภายในโครงการเป็น Python source code ไม่ต้องติดตั้งเป็นโปรแกรมแยก ส่วน
NOVOPlasty, GetOrganelle, MIA, Pilon และ MITObim ดาวน์โหลดเป็นแพ็กเกจจาก
[Bioconda](https://bioconda.github.io/) โดยใช้ dependencies จาก Bioconda และ
conda-forge ผ่าน `micromamba`

## 2. ข้อมูลที่ใช้: แหล่งที่มา การดาวน์โหลด และตำแหน่งจัดเก็บ

ข้อมูลไม่ใช่ software จึงไม่ได้ “ติดตั้ง” ลงใน micromamba environment แต่ดาวน์โหลด
มาเก็บใต้ `asem-refine-benchmarks/data/` ซึ่งถูกแยกออกจาก source code และไม่ควร
สมมติว่าจะมีอยู่หลัง clone repository ใหม่ ต้องรันสคริปต์เตรียมข้อมูลก่อน

### 2.1 Reference genomes จาก NCBI

`scripts/fetch_data.py` ดาวน์โหลด nucleotide FASTA จากฐานข้อมูล NCBI
Nucleotide ผ่าน NCBI E-utilities โดยใช้ accession ที่กำหนดไว้ในสคริปต์:

```bash
cd /home/punpiti/OneDrive/genome/asem-refine-benchmarks/scripts
micromamba run -n genome python fetch_data.py
```

endpoint ที่ใช้คือ:

```text
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
```

request ใช้ `db=nuccore`, `rettype=fasta` และ `retmode=text` สคริปต์พัก 0.34
วินาทีระหว่าง request เพื่อไม่เกินอัตรา unauthenticated NCBI E-utilities ประมาณ
3 requests ต่อวินาที ถ้าไฟล์ปลายทางมีอยู่และไม่ว่าง สคริปต์จะข้ามไฟล์นั้น จึงรัน
ซ้ำได้โดยไม่ดาวน์โหลดใหม่ทั้งหมด

ตำแหน่งจัดเก็บ default:

```text
data/phase1_reproduction/       mtDNA 6 species สำหรับ main simulation grid
data/ieee2021_mtdna_panel/      mtDNA 30 species สำหรับ divergence panel
data/refs/                      references สำหรับ real-WGS validation
```

ชุดหลัก 6 species ได้แก่:

| Species | NCBI accession | บทบาท |
|---|---|---|
| *Homo sapiens* | NC_012920 | Phase-1 และ human rCRS evaluation target |
| *Gorilla gorilla* | NC_001645 | Phase-1 same-family comparison |
| *Saimiri boliviensis* | NC_021966 | Phase-1 same-genus comparison |
| *Saimiri sciureus* | NC_012775 | Phase-1 same-genus comparison |
| *Aotus azarai* | NC_021939 | Phase-1 primate reference/target |
| *Varecia variegata* | NC_012773 | Phase-1 และ divergent real-read starting reference |

ชุด 30 species ใช้ accessions ที่แสดงใน
`docs/ieee2021_mtdna_panel_manifest.md` และใน `MANIFESTS` ของ
`scripts/fetch_data.py` ไฟล์ที่ใช้อยู่จริงมีขนาดประมาณ 16.6–17.3 kb ต่อ genome

reference เพิ่มเติมสำหรับ real-WGS ได้แก่:

| ไฟล์ | Accession/ที่มา | ใช้เป็น |
|---|---|---|
| `data/refs/chrM_rCRS.fasta` | Human NC_012920 | evaluation truth/target |
| `data/refs/chimp_NC_001643.fasta` | Chimpanzee NC_001643 | same-family starting reference |
| `data/refs/chrM_for_cram.fasta` | สำเนา NC_012920 ที่เปลี่ยน header เป็น `chrM` | reference สำหรับอ่าน remote CRAM |

### 2.2 การจำลอง reads สำหรับ benchmark หลัก

simulation ไม่ได้ดาวน์โหลด read files สำเร็จรูป แต่สร้าง reads ใหม่จาก target
mtDNA ในแต่ละ job ด้วย `simulate_shotgun_reads()`:

- single-end read length 150 bp
- coverage 1–8X สำหรับ main ASEM/SRSC/Hybrid grid
- NOVOPlasty/GetOrganelle มี depth เพิ่มเป็น 10, 15, 20 และ 30X
- 3 replicates ใน full grid
- ใช้ `stable_seed(target_name, depth, replicate)`
- methods ที่มี species/depth/replicate เดียวกันจึงได้ reads เหมือนกันแบบ
  read-for-read

main grid สร้างทุก ordered pair ที่ reference ไม่เท่ากับ target จาก 6 species:
30 pairs × 8 depths × 3 replicates = 720 jobs ต่อวิธี

reads ถูกสร้างใน memory สำหรับวิธีภายใน ส่วน external wrappers เขียนเป็น FASTA
หรือ FASTQ ชั่วคราวใน `/tmp` แยกต่อ job และลบเมื่อ job จบ จึงไม่มี simulated-read
dataset ขนาดใหญ่ติดตั้งค้างไว้ใน `data/`

### 2.3 Whole-genome sequencing source สำหรับ real-WGS validation

real-read validation เริ่มจากข้อมูล whole-genome sequencing ของมนุษย์ แต่
benchmark **ไม่ได้ประกอบหรือปรับปรุง nuclear whole genome** ข้อมูลต้นทางคือ
ตัวอย่าง NA07000 จาก 1000 Genomes Project 30× resequencing collection:

```text
ENA project:  PRJEB31736
ENA run:      ERR3239279
Sample:       NA07000
Source file:  NA07000.final.cram
Remote URL:   https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239279/NA07000.final.cram
Full size:    ประมาณ 18 GB
```

เพื่อหลีกเลี่ยงการดาวน์โหลด CRAM ทั้ง genome สคริปต์ใช้ remote indexed range
query ของ samtools ดึงเฉพาะ records ที่ map กับ `chrM`:

```bash
cd /home/punpiti/OneDrive/genome/asem-refine-benchmarks

# สร้าง reference ที่ header ตรงกับ @SQ SN:chrM ใน CRAM
awk 'NR==1 {print ">chrM"; next} {print}' \
  data/refs/chrM_rCRS.fasta > data/refs/chrM_for_cram.fasta
micromamba run -n genome samtools faidx data/refs/chrM_for_cram.fasta

# ดาวน์โหลด CRAI และดึงเฉพาะ chrM slice เป็น BAM
cd scripts
SAMTOOLS_BIN=/home/punpiti/.local/share/mamba/envs/genome/bin/samtools \
./extract_real_wgs_reads.sh
```

คำสั่งหลักภายในสคริปต์คือ:

```bash
curl -sf "${CRAM_URL}.crai" -o data/NA07000.final.cram.crai

samtools view -b \
  --reference data/refs/chrM_for_cram.fasta \
  -X "$CRAM_URL" data/NA07000.final.cram.crai \
  chrM -o data/NA07000_chrM.bam
```

ไฟล์ที่ได้และตำแหน่งจริงในโครงการ:

```text
data/NA07000.final.cram.crai       1,526,808 bytes
data/NA07000_chrM.bam             91,509,533 bytes
data/refs/chrM_for_cram.fasta
data/refs/chrM_for_cram.fasta.fai
```

การรันที่รายงานได้ 2,254,405 records ใน chrM BAM ก่อน filtering จากนั้น
`run_real_wgs*.py` เรียก:

```bash
samtools view -F 3332 data/NA07000_chrM.bam
```

flag mask 3332 ตัด unmapped (`0x4`), secondary (`0x100`), duplicate (`0x400`)
และ supplementary (`0x800`) records แล้วเก็บเฉพาะ sequence ที่ยาว exactly
150 bp เหลือ usable pool 1,860,886 reads

หมายเหตุเรื่อง strand: SEQ ของ mapped reverse-strand records ใน BAM อยู่ใน
reference orientation อยู่แล้ว runner จึงไม่ reverse-complement ซ้ำ

### 2.4 การเก็บข้อมูลไว้นอก repository

default คือ `asem-refine-benchmarks/data/` แต่หากต้องเก็บ BAM บน disk อื่น ให้ส่ง
output directory ให้ extraction script และชี้ runner ด้วย `GENOME_WGS_DATA_DIR`:

```bash
cd /home/punpiti/OneDrive/genome/asem-refine-benchmarks/scripts

CHRM_REF_FASTA=/path/to/data/refs/chrM_for_cram.fasta \
SAMTOOLS_BIN=/path/to/samtools \
./extract_real_wgs_reads.sh /path/to/data

GENOME_WGS_DATA_DIR=/path/to/data \
SAMTOOLS_BIN=/path/to/samtools \
micromamba run -n genome python run_real_wgs_benchmark.py
```

directory ภายนอกต้องรักษา layout นี้:

```text
/path/to/data/NA07000_chrM.bam
/path/to/data/refs/chrM_rCRS.fasta
/path/to/data/refs/chimp_NC_001643.fasta
```

สำหรับ hybrid runner ไฟล์ Varecia ยังอ่านจาก
`data/phase1_reproduction/Varecia_variegata_NC_012773.fasta` ภายใน repository

## 3. Environment หลัก

repository มี `setup_env.sh` สำหรับสร้าง environment หลักชื่อ `asem-bench`:

```bash
cd /home/punpiti/OneDrive/genome/asem-refine-benchmarks
bash setup_env.sh
```

คำสั่งหลักภายในคือ:

```bash
micromamba create -n asem-bench -y \
  -c bioconda -c conda-forge \
  --channel-priority flexible \
  "samtools>=1.10" scikit-bio numpy pandas matplotlib seaborn scipy

micromamba run -n asem-bench pip install parasail
```

การทดลองล่าสุดในเครื่องนี้เรียก Python จาก environment `genome` ซึ่งมี Python
3.14.6, samtools 1.21, scikit-bio 0.7.3, parasail 1.3.4 และ NOVOPlasty 4.3.5
คำสั่งต่อไปจึงใช้ `micromamba run -n genome` เพื่อแสดงการรันจริง หากสร้างเครื่อง
ใหม่ด้วย `setup_env.sh` ให้เปลี่ยนเป็น `-n asem-bench` และติดตั้ง NOVOPlasty เพิ่ม

## 4. วิธีภายในโครงการ

### ASEM non-recursive

```bash
cd /home/punpiti/OneDrive/genome/asem-refine-benchmarks
ASEM_GRID_WORKERS=14 \
micromamba run -n genome python scripts/run_grid_baseline1.py
```

### ASEM recursive

```bash
cd /home/punpiti/OneDrive/genome/asem-refine-benchmarks
ASEM_GRID_WORKERS=14 \
micromamba run -n genome python scripts/run_grid_baseline2.py
```

### SRSC

```bash
cd /home/punpiti/OneDrive/genome/asem-refine-benchmarks
micromamba run -n genome python scripts/run_grid_baseline3.py
```

สามวิธีนี้ใช้ coverage 1–8X, 3 replicates และ 30 ordered pairs รวม 720 jobs
ต่อวิธี `ASEM_GRID_WORKERS` กำหนดจำนวน jobs ที่ทำพร้อมกัน; default 14 สำหรับ
ASEM runners ส่วน ASEM recursive วน refinement ได้ไม่เกิน 6 iterations

### ASEM-Hybrid

runner รองรับทั้ง `legacy` และ `boundary` แต่ default ใน source ปัจจุบันยังเป็น
`legacy` ขณะที่ manuscript ล่าสุดใช้ boundary-recruited variant จึงต้องระบุ:

```bash
cd /home/punpiti/OneDrive/genome/asem-refine-benchmarks
ASEM_HYBRID_VARIANT=boundary \
ASEM_HYBRID_WORKERS=14 \
micromamba run -n genome python scripts/run_grid_baseline_hybrid.py
```

หากไม่กำหนด `ASEM_HYBRID_VARIANT=boundary` ผลจะไม่ตรงกับ manuscript ล่าสุด

### Full rerun หลังแก้ tie-breaking (revision 2026-09-30)

การยืนยันผลหลังเปลี่ยนกฎจาก alphabetic `argmax()` เป็นการคงฐานเดิมเมื่อ
nucleotide เสมอกัน และไม่รับ tied candidate insertion ใช้ shell script เดียว:

```bash
cd /home/punpiti/OneDrive/genome/asem-refine-benchmarks
bash scripts/run_post_tie_rule_grids.sh
```

สคริปต์ทำงานตามลำดับนี้โดยอัตโนมัติ:

1. รัน unit tests ทั้งหมด
2. รัน ASEM non-recursive 720 jobs
3. รัน ASEM recursive 720 jobs
4. รัน boundary ASEM-Hybrid 720 jobs
5. ตรวจความครบ 720 jobs ต่อวิธีและสร้างผลเปรียบเทียบก่อน/หลังแก้กฎ

default ใช้ 14 workers หากต้องการเปลี่ยนจำนวน core:

```bash
ASEM_GRID_WORKERS=10 bash scripts/run_post_tie_rule_grids.sh
```

runner เขียนผลเมื่อแต่ละ job จบและอ่าน job keys เดิมก่อนเริ่ม จึง **resume ได้**:
หาก process ถูกหยุด ให้เรียกคำสั่งเดิมซ้ำโดยไม่ต้องลบ CSV และจะข้าม jobs ที่เสร็จแล้ว
ห้ามลบหรือเปลี่ยนชื่อไฟล์ `grid_results.post_tie_rule_partial_20260929.csv`
ระหว่างการ resume

หากต้องการปล่อยให้รันโดยไม่เปิด terminal ค้าง:

```bash
cd /home/punpiti/OneDrive/genome/asem-refine-benchmarks
mkdir -p scripts/results/post_tie_rule_logs
nohup bash scripts/run_post_tie_rule_grids.sh \
  > scripts/results/post_tie_rule_logs/controller.log 2>&1 &
echo $! > scripts/results/post_tie_rule_logs/controller.pid
```

ตรวจสถานะได้โดยไม่ต้องให้ AI อ่าน log:

```bash
bash scripts/check_post_tie_rule_status.sh
```

logs แยกตามวิธีอยู่ใน:

```text
scripts/results/post_tie_rule_logs/baseline1.log
scripts/results/post_tie_rule_logs/baseline2.log
scripts/results/post_tie_rule_logs/hybrid_boundary.log
```

ผล post-tie ที่กำลังสร้างอยู่ใน:

```text
scripts/results/baseline1_ieee_access/grid_results.post_tie_rule_partial_20260929.csv
scripts/results/baseline2_ojemb/grid_results.post_tie_rule_partial_20260929.csv
scripts/results/baseline_hybrid_boundary/grid_results.post_tie_rule_partial_20260929.csv
```

เมื่อครบ สคริปต์ `summarize_post_tie_rule.py` จะสร้าง:

```text
scripts/results/post_tie_rule_summary.json
```

ไฟล์นี้รายงานจำนวน jobs/rows, งานที่ผลสุดท้ายเปลี่ยน, delta ของ F1/identity/
recall/precision และค่าที่ต้องใช้ตรวจตัวเลขใน manuscript การรัน grid ทั้งหมดเป็น
งานของ shell/Python scripts; AI จำเป็นเฉพาะตอนอ่าน summary แล้วแก้ข้อความและตาราง
ใน manuscript เท่านั้น

## 5. NOVOPlasty

### ดาวน์โหลดและติดตั้ง

ติดตั้งแพ็กเกจ `novoplasty` จาก Bioconda:

```bash
micromamba install -n genome -y \
  -c bioconda -c conda-forge \
  novoplasty=4.3.5
```

แหล่งแพ็กเกจ: <https://anaconda.org/bioconda/novoplasty>

executable ที่ใช้จริง:

```text
/home/punpiti/.local/share/mamba/envs/genome/bin/NOVOPlasty4.3.5.pl
```

### การเรียกและ parameter

`scripts/ext_novoplasty.py` สร้าง scratch directory ต่อ job เขียน reads เป็น
FASTQ และเรียก:

```bash
perl /home/punpiti/.local/share/mamba/envs/genome/bin/NOVOPlasty4.3.5.pl \
  -c config.txt
```

config ที่ใช้:

```text
Type                    = mito
Genome Range            = 12000-22000
K-mer                   = 33
Max memory              = 4 GB
Seed Input              = seed_ref.fasta
Reference sequence      = seed_ref.fasta
Read Length             = 150
Insert size             = 300
Platform                = illumina
Single/Paired           = SE
Combined reads          = reads.fastq
Insert size auto        = yes
Use Quality Scores      = no
```

wrapper เลือก sequence ที่ยาวที่สุดจาก Contigs, Circularized, Merged หรือ Option
FASTA เป็น output assembly

### การรัน grid

```bash
cd /home/punpiti/OneDrive/genome/asem-refine-benchmarks
micromamba run -n genome python scripts/run_grid_ext_novoplasty.py
```

- coverage 1–8, 10, 15, 20 และ 30X
- 3 replicates ต่อ ordered pair/depth
- timeout 300 วินาทีต่อ job
- รันพร้อมกันสูงสุด 14 jobs ตามค่าที่กำหนดใน runner

NOVOPlasty เป็น de novo organelle assembler ผลล้มเหลวที่ 1–8X จึงต้องอธิบายว่า
เป็นผลภายใต้ ultra-low coverage และ configuration นี้ ไม่ใช่ข้อสรุปทั่วไปว่า
NOVOPlasty ใช้งานไม่ได้

## 6. GetOrganelle

### ดาวน์โหลดและติดตั้ง

```bash
micromamba create -n getorganelle -y \
  -c bioconda -c conda-forge \
  getorganelle=1.7.7.0
```

แหล่งแพ็กเกจ: <https://anaconda.org/bioconda/getorganelle>

environment ที่ใช้จริงมี GetOrganelle 1.7.7.0, SPAdes 3.15.5 และ Bowtie2 2.5.4
ตัวโปรแกรมอยู่ที่:

```text
/home/punpiti/.local/share/mamba/envs/getorganelle/bin/get_organelle_from_reads.py
```

### การเรียกและ parameter

`scripts/ext_getorganelle.py` เพิ่ม environment นี้เข้า `PATH` เพื่อให้หา Bowtie2
และ SPAdes ที่ตรงกัน แล้วเรียก:

```bash
/home/punpiti/.local/share/mamba/envs/getorganelle/bin/python \
  /home/punpiti/.local/share/mamba/envs/getorganelle/bin/get_organelle_from_reads.py \
  -u reads.fastq \
  -s seed_ref.fasta \
  -F animal_mt \
  -o out \
  -R 10 \
  -k 21,45,65,85,105 \
  -t 1
```

- `-u`: single-end reads
- `-s`: starting seed
- `-F animal_mt`: animal mitochondrial genome
- `-R 10`: 10 extension rounds
- `-k 21,45,65,85,105`: ชุด k-mer
- `-t 1`: หนึ่ง thread ต่อ job เพราะ grid รันหลาย jobs พร้อมกัน

wrapper เลือก sequence ที่ยาวที่สุดจาก `*path_sequence.fasta`

### การรัน grid

```bash
cd /home/punpiti/OneDrive/genome/asem-refine-benchmarks
GETORGANELLE_ENV_BIN=/home/punpiti/.local/share/mamba/envs/getorganelle/bin \
micromamba run -n genome python scripts/run_grid_ext_getorganelle.py
```

ใช้ coverage 1–8, 10, 15, 20 และ 30X, 3 replicates, timeout 600 วินาทีต่อ job
และอนุญาตสูงสุด 14 jobs พร้อมกัน

## 7. MIA (Mapping Iterative Assembler)

### ดาวน์โหลดและติดตั้ง

```bash
micromamba create -n mia-assembler -y \
  -c conda-forge -c bioconda \
  --channel-priority flexible \
  mapping-iterative-assembler=1.0
```

แหล่งแพ็กเกจ: <https://anaconda.org/bioconda/mapping-iterative-assembler>

ไฟล์แพ็กเกจจริงมาจาก `conda.anaconda.org/bioconda/linux-64/` ชื่อ
`mapping-iterative-assembler-1.0-h503566f_7.conda`

executables:

```text
/home/punpiti/.local/share/mamba/envs/mia-assembler/bin/mia
/home/punpiti/.local/share/mamba/envs/mia-assembler/bin/ma
```

### การเรียกและ parameter

`scripts/ext_mia.py` เขียน reference และ reads เป็น FASTA แล้วเรียก:

```bash
mia \
  -r reference.fasta \
  -f reads.fasta \
  -m assembly \
  -c -D -F
```

- `-r`: starting reference
- `-f`: input reads
- `-m`: output root
- `-c`: circular-reference mode
- `-D`: distant-reference mode; คง low-scoring reads ระหว่าง iterative assembly
- `-F`: เขียนเฉพาะ final MALN assembly

จากนั้นเลือก `assembly.<iteration>` หมายเลขสูงสุดและแปลง MALN เป็น FASTA:

```bash
ma -M assembly.<final_iteration> -f 5 -I MIA_consensus
```

ข้อควรระวัง: `ma -f 5` หมายถึง output format 5 หรือ FASTA ไม่ใช่ minimum
mapping quality 5 ถ้า manuscript ใช้คำอธิบายหลังต้องแก้ให้ตรงกับคำสั่งจริง

### การรัน representative benchmark

```bash
cd /home/punpiti/OneDrive/genome/asem-refine-benchmarks
MIA_ENV_BIN=/home/punpiti/.local/share/mamba/envs/mia-assembler/bin \
MIA_SCOPE=representative \
MIA_TIMEOUT_S=300 \
MIA_PARALLEL_RUNS=6 \
micromamba run -n genome python scripts/run_grid_ext_mia.py
```

representative subset มี 12 jobs:

```text
Saimiri_boliviensis -> Saimiri_sciureus    same genus
Homo_sapiens        -> Gorilla_gorilla      same family
Homo_sapiens        -> Varecia_variegata    same order
coverage: 1, 2, 4 และ 8X
replicate: 0
```

ต้องกำหนด `MIA_TIMEOUT_S=300` เพราะผลที่เก็บไว้และ manuscript ใช้ cutoff 300
วินาที แม้ default ใน runner ปัจจุบันยังเป็น 180 วินาที

full 720-job grid เรียกได้ด้วย:

```bash
MIA_SCOPE=full \
MIA_TIMEOUT_S=300 \
MIA_PARALLEL_RUNS=6 \
micromamba run -n genome python scripts/run_grid_ext_mia.py
```

แต่ MIA ใช้เวลามากหรือ timeout ได้แม้กับ mitochondrial genome จึงใช้ matched
representative subset ในการตอบ reviewer

## 8. Pilon

### ดาวน์โหลดและติดตั้ง

```bash
micromamba create -n pilon-benchmark -y \
  -c conda-forge -c bioconda \
  pilon=1.24 bwa=0.7.19 samtools=1.24
```

แหล่งแพ็กเกจ:

- <https://anaconda.org/bioconda/pilon>
- <https://anaconda.org/bioconda/bwa>
- <https://anaconda.org/bioconda/samtools>

ตำแหน่งที่ wrapper ใช้:

```text
/home/punpiti/.local/share/mamba/envs/pilon-benchmark/bin/bwa
/home/punpiti/.local/share/mamba/envs/pilon-benchmark/bin/samtools
/home/punpiti/.local/share/mamba/envs/pilon-benchmark/bin/java
/home/punpiti/.local/share/mamba/envs/pilon-benchmark/share/pilon-1.24-0/pilon.jar
```

### การเรียกและ parameter

Pilon รับ BAM ที่ map กับ starting reference ดังนั้น `scripts/ext_pilon.py` รัน:

```bash
bwa index reference.fasta
bwa mem -t 1 reference.fasta reads.fasta > aligned.sam
samtools sort -o aligned.bam aligned.sam
samtools index aligned.bam

java -Xms512m -Xmx1g \
  -jar pilon.jar \
  --genome reference.fasta \
  --unpaired aligned.bam \
  --output pilon \
  --outdir output_directory \
  --fix all
```

- `bwa mem -t 1`: หนึ่ง mapping thread ต่อ job
- `--unpaired`: reads จำลองเป็น single-end
- `--fix all`: เปิด correction categories ที่ Pilon รองรับทั้งหมด
- `-Xms512m -Xmx1g`: Java heap 512 MB–1 GB ต่อ job

### การรัน representative benchmark

```bash
cd /home/punpiti/OneDrive/genome/asem-refine-benchmarks
PILON_PARALLEL_RUNS=6 \
PILON_TIMEOUT_S=180 \
micromamba run -n genome python scripts/run_grid_ext_pilon.py
```

ใช้สาม representative pairs, coverage 1, 2, 4 และ 8X, replicate 0 เช่นเดียวกับ
MIA รวม 12 jobs

## 9. MITObim: ติดตั้งแล้วแต่ไม่รวมใน benchmark

```bash
micromamba create -n mitobim -y \
  -c bioconda -c conda-forge \
  mitobim=1.9.1
```

แหล่งแพ็กเกจ: <https://anaconda.org/bioconda/mitobim>

environment มี MITObim 1.9.1, MIRA 4.0.2, Perl 5.22 และ Python 2.7
executable อยู่ที่:

```text
/home/punpiti/.local/share/mamba/envs/mitobim/bin/MITObim.pl
```

แต่ `mirabait` ที่ bundle มากับ MIRA ใช้ legacy `vsyscall` และเกิด segmentation
fault บน WSL2 kernel 6.6 ที่ใช้ ปัญหาเกิดแม้กับ one-read smoke test และทดสอบ
MIRA 4.9.6 แยกแล้ว จึงไม่มี MITObim grid result นี่เป็น technical exclusion
ไม่ใช่ผลว่า MITObim ล้มเหลวเพราะ coverage ต่ำ

## 10. Real-WGS benchmark และการทดลองประกอบอื่น

### 10.1 Original real-WGS benchmark

การทดลองนี้ใช้ chimpanzee mtDNA `NC_001643` เป็น starting reference, human rCRS
`NC_012920` เป็น evaluation target และ subsample reads จาก filtered NA07000 pool
โดยไม่คืน reads ที่สุ่มไปแล้ว:

```bash
cd /home/punpiti/OneDrive/genome/asem-refine-benchmarks
GENOME_WGS_DATA_DIR=/home/punpiti/OneDrive/genome/asem-refine-benchmarks/data \
SAMTOOLS_BIN=/home/punpiti/.local/share/mamba/envs/genome/bin/samtools \
micromamba run -n genome python scripts/run_real_wgs_benchmark.py
```

parameter:

- methods: ASEM non-recursive และ ASEM recursive
- starting reference: chimpanzee mtDNA NC_001643
- truth: human rCRS NC_012920
- depths: 1–8X
- 3 replicates
- read length: exactly 150 bp
- maximum iterations: 6
- จำนวน reads ต่อ job: `round(depth × 16,569 / 150)` ได้แก่ประมาณ 110 reads
  ที่ 1X ถึง 884 reads ที่ 8X

ผลเก็บที่:

```text
scripts/results/real_wgs_benchmark/grid_results.csv
```

### 10.2 Boundary-Hybrid real-WGS benchmark

การทดลองล่าสุดใช้ read subsets เดียวกันเปรียบเทียบ starting references สองระดับ:

1. chimpanzee NC_001643 → human: same-family, ใช้ตรวจว่า Hybrid ไม่ทำลายกรณีที่
   ordinary mapping ทำงานได้ดีอยู่แล้ว
2. *Varecia variegata* NC_012773 → human: same-order, ใช้ทดสอบประโยชน์ของ
   boundary recruitment เมื่อ reference ห่างกว่า

คำสั่งที่ตรงกับ manuscript ล่าสุด:

```bash
cd /home/punpiti/OneDrive/genome/asem-refine-benchmarks
GENOME_WGS_DATA_DIR=/home/punpiti/OneDrive/genome/asem-refine-benchmarks/data \
SAMTOOLS_BIN=/home/punpiti/.local/share/mamba/envs/genome/bin/samtools \
ASEM_HYBRID_VARIANT=boundary \
ASEM_REAL_WGS_WORKERS=2 \
micromamba run -n genome python scripts/run_real_wgs_hybrid_benchmark.py
```

parameter เพิ่มเติม:

- ASEM mapping threshold `tau = 0.5`
- `boundary_flank = 150` สำหรับ boundary variant
- maximum iterations 6
- 1 worker ภายในแต่ละ algorithm job
- `ASEM_REAL_WGS_WORKERS=2` ควบคุมจำนวน jobs พร้อมกัน
- Pan reference รัน ASEM non-recursive, ASEM recursive และ ASEM-Hybrid
- Varecia reference รัน ASEM non-recursive และ ASEM-Hybrid

ผลเก็บที่:

```text
scripts/results/real_wgs_hybrid_boundary/final_results.csv
```

ห้ามละ `ASEM_HYBRID_VARIANT=boundary` เพราะ default ของ runner ยังเป็น `legacy`
และจะเขียนไปคนละ results directory

### 10.3 ขอบเขตของคำว่า whole/full genome

ข้อมูล NA07000 ต้นทางเป็น whole-genome sequencing จริง แต่ benchmark นี้ใช้เพียง
mitochondrial-mapped read slice และประเมิน mitochondrial reference ยาวประมาณ
16.6 kb เท่านั้น ไม่ได้ทดสอบการ reconstruct หรือ refine nuclear genome ขนาด
ประมาณ 3 Gb ดังนั้นใน manuscript ควรเรียกการทดลองนี้ว่า “real-WGS-derived
mitochondrial-read validation” หรือ “real-read mtDNA validation” และไม่ควรเรียก
ว่า “full-genome benchmark”

ถ้าต้องการอ้างความสามารถระดับ whole nuclear genome จะต้องออกแบบ benchmark ใหม่
แยกต่างหาก รวมถึง chromosome/region selection, memory scaling, repeat handling,
structural variation และ evaluation truth ซึ่งยังไม่มีใน repository นี้

### 10.4 การทดลองประกอบอื่น

```bash
cd /home/punpiti/OneDrive/genome/asem-refine-benchmarks

# Read-length experiment
micromamba run -n genome python scripts/run_read_length_experiment.py

# Computational-cost profiling
micromamba run -n genome python scripts/profile_computational_cost.py

# Reviewer-requested diagnostics
micromamba run -n genome python scripts/run_parameter_sensitivity.py
micromamba run -n genome python scripts/run_circular_boundary_diagnostic.py
micromamba run -n genome python scripts/run_evaluator_endgap_sanity.py
```

## 11. ผลลัพธ์และการ resume

runner เขียน CSV ใต้:

```text
scripts/results/<benchmark-name>/grid_results.csv
```

ตัวอย่าง:

```text
scripts/results/baseline1_ieee_access/grid_results.csv
scripts/results/baseline2_ojemb/grid_results.csv
scripts/results/baseline_hybrid_boundary/grid_results.csv
scripts/results/ext_mia_representative/grid_results.csv
scripts/results/ext_pilon_representative/grid_results.csv
```

หลาย runner รองรับการ resume โดยอ่าน CSV เดิมและข้าม `(reference, target,
depth, replicate)` ที่มีอยู่แล้ว การเรียกซ้ำจึงอาจไม่มี pending jobs และไม่ใช่
fresh rerun หากต้องการรันใหม่ทั้งหมดควรทำในสำเนา repository หรือย้าย results เดิม
ไปเก็บก่อน ไม่ควรเขียนทับ raw results ของ manuscript โดยไม่มีสำเนา

สคริปต์ `make_*.py` อ่าน raw CSV เพื่อสร้าง figures/tables เช่น:

```bash
micromamba run -n genome python scripts/make_external_tools_depth_chart_english.py
micromamba run -n genome python scripts/make_hybrid_comparison_figure.py
micromamba run -n genome python scripts/make_computational_cost_table.py
micromamba run -n genome python scripts/make_report.py
```

## 12. สถานะ reproducibility ของ revision snapshot

revision snapshot `v1.2.0` รวม tie-safe consensus, full post-tie grids,
MIA/Pilon wrappers และ raw results, parameter sensitivity, circular-boundary
diagnostic, evaluator sanity check และ strand/origin-normalized external-tool
evaluator แล้ว คำสั่ง Hybrid ระบุ `ASEM_HYBRID_VARIANT=boundary` และเอกสาร MIA
อธิบาย `ma -f 5` ว่าเป็น FASTA output format พร้อมระบุ timeout 300 วินาที

ขั้น release ภายนอกที่ต้องทำหลัง commit คือสร้าง tag/GitHub release และ Zenodo
version ใหม่ จากนั้นจึงแทน version, commit hash และ archived DOI ใน Code
availability ของ revised manuscript ห้ามอ้าง v1.1.0/commit เดิมว่าเป็น snapshot
ของผล revision

## 13. ตรวจสอบการติดตั้ง

```bash
micromamba env list

/home/punpiti/.local/share/mamba/envs/genome/bin/NOVOPlasty4.3.5.pl --help
/home/punpiti/.local/share/mamba/envs/getorganelle/bin/get_organelle_from_reads.py --version
/home/punpiti/.local/share/mamba/envs/mia-assembler/bin/mia
/home/punpiti/.local/share/mamba/envs/mia-assembler/bin/ma
/home/punpiti/.local/share/mamba/envs/pilon-benchmark/bin/java \
  -jar /home/punpiti/.local/share/mamba/envs/pilon-benchmark/share/pilon-1.24-0/pilon.jar \
  --version
```

โปรแกรมรุ่นเก่าบางตัวอาจไม่มี `--version` หรือคืน exit status ไม่เป็นศูนย์เมื่อ
เรียก help จึงควรตรวจทั้ง executable output และ package metadata ใต้
`<environment>/conda-meta/`
