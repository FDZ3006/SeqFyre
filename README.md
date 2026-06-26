# SeqFyre

Pipeline analisis GC content untuk sekuens nukleotida. Seluruh logika inti diimplementasikan dalam Python. Membaca berkas FASTA atau FASTQ, menghitung frekuensi nukleotida, mengurutkan sekuens berdasarkan GC content, dan menghasilkan visualisasi serta ekspor CSV. HTML, CSS, dan JavaScript yang ada di repositori hanya berfungsi sebagai antarmuka tampilan web dan tidak mengandung logika analisis.

**Akses langsung (tanpa instalasi):** https://web-production-4b2df.up.railway.app

**Download aplikasi desktop:** [SeqFyre.exe - GitHub Releases](https://github.com/FDZ3006/seqfyre/releases/latest) (Windows, tidak perlu install Python)

---

**Mata Kuliah:** Struktur Data Bioinformatika (BIF1222)  
**Program Studi:** Bioinformatika, IPB University  
**Penulis:** Fajar Dhuha Zafran (G0401241009)

---

## Daftar Isi

- [Latar Belakang](#latar-belakang)
- [Dataset](#dataset)
- [Fitur](#fitur)
- [Struktur Proyek](#struktur-proyek)
- [Cara Menggunakan](#cara-menggunakan)
  - [Via Website](#via-website-disarankan)
  - [Via Desktop App](#via-desktop-app)
  - [Via Lokal](#via-lokal)
  - [Via Command Line](#via-command-line)
- [Hasil Analisis](#hasil-analisis)
- [Eksplorasi Data](#eksplorasi-data)
- [Pengujian](#pengujian)
- [English Documentation](#english-documentation)

---

## Latar Belakang

Stabilitas struktural asam nukleat dipengaruhi oleh komposisi basa penyusunnya. Pasangan basa Guanin-Sitosin (G-C) memiliki tiga ikatan hidrogen dan gaya base-stacking yang lebih kuat dibandingkan pasangan Adenin-Timin (A-T) yang hanya memiliki dua ikatan hidrogen. Perbedaan ini menyebabkan molekul dengan kandungan GC tinggi lebih tahan terhadap suhu tinggi.

Bakteri termofilik, yang hidup di lingkungan bersuhu ekstrem seperti sumber air panas dan ventilasi hidrotermal, memanfaatkan mekanisme ini secara evolusioner dengan memperkaya kandungan G dan C pada gen struktural RNA-nya. Sebaliknya, bakteri mesofilik yang hidup pada suhu sedang memiliki kandungan GC yang lebih rendah.

Gen 16S rRNA dipilih sebagai dataset karena merupakan biomarker filogenetik universal dengan panjang sekitar 1.500 bp, cukup terkonservasi untuk perbandingan lintas spesies namun tetap mencerminkan tekanan seleksi suhu pada komposisi nukleotidanya.

Ketika pipeline ini mengurutkan sekuens berdasarkan GC content secara menurun, bakteri termofilik secara konsisten muncul di posisi teratas. Hal ini memvalidasi bahwa pipeline tidak sekadar melakukan operasi komputasi, tetapi juga mencerminkan fenomena biologis yang nyata.

---

## Dataset

Berkas: `data/dataset_16S_rRNA_SeqFyre.fasta`

10 sekuens gen 16S rRNA dari NCBI RefSeq, terdiri dari 4 bakteri termofilik dan 6 bakteri mesofilik.

| Aksesi | Organisme | Kelompok | Suhu Optimal |
|--------|-----------|----------|-------------|
| NR_075056.2 | Aquifex aeolicus VF5 | Termofilik | ~95°C |
| NR_037066.1 | Thermus thermophilus HB8 | Termofilik | ~70-75°C |
| NR_102775.2 | Thermotoga maritima MSB8 | Termofilik | ~80°C |
| NR_115284.2 | Geobacillus stearothermophilus | Termofilik | ~65°C |
| NR_112116.2 | Bacillus subtilis | Mesofilik | ~37°C |
| NR_114042.1 | Escherichia coli | Mesofilik | ~37°C |
| NR_024570.1 | Escherichia coli (strain U 5/41) | Mesofilik | ~37°C |
| NR_026078.1 | Pseudomonas aeruginosa | Mesofilik | ~37°C |
| NR_112081.1 | Herbaspirillum rubrisubalbicans | Mesofilik | ~30°C |
| NR_118997.2 | Staphylococcus aureus | Mesofilik | ~37°C |

---

## Fitur

Seluruh poin berikut diimplementasikan dalam Python menggunakan struktur data List dan Dictionary sebagai inti pipeline.
Poin wajib yang dipenuhi:

- Membaca berkas FASTA atau FASTQ
- Menyimpan data dalam struktur List
- Menghitung frekuensi nukleotida menggunakan Dictionary
- Mengurutkan sekuens berdasarkan GC content
- Menampilkan 3 sekuens terbaik (dapat dikonfigurasi: 3, 5, atau 10)
- Visualisasi grafik berdasarkan nilai GC
- Menulis hasil ke berkas CSV

Fitur tambahan:

- Antarmuka web (Flask) yang dapat diakses langsung dari browser tanpa instalasi
- Dukungan unggah berkas FASTA, FASTQ, dan ZIP
- 4 grafik eksplorasi data (histogram, bar chart, komposisi nukleotida, scatter plot)
- Unduhan hasil analisis dalam format ZIP berisi CSV dan grafik
- Kode berorientasi objek dengan kelas terpisah (SequenceRecord, Parser, Analyzer)
- Antarmuka dwibahasa (Indonesia dan Inggris)
- Tiga pilihan tema tampilan (terang, gelap, Fyre)

---

## Struktur Proyek

```
seqfyre/
├── seqfyre/               # Paket inti (OOP)
│   ├── models.py          # SequenceRecord - model satu sekuens
│   ├── parser.py          # Parser - membaca FASTA/FASTQ/ZIP
│   ├── analyzer.py        # Analyzer - pengurutan, statistik, grafik, CSV
│   └── packager.py        # Pengemas hasil menjadi ZIP
├── app.py                 # Aplikasi web Flask
├── cli.py                 # Pipeline berbasis command line
├── templates/index.html   # Antarmuka web
├── static/                # CSS dan JavaScript
├── data/                  # Dataset 16S rRNA
├── tests/                 # Pengujian unit
├── assets/                # Grafik untuk dokumentasi
└── requirements.txt
```

### Penjelasan Kelas

| Kelas | Tanggung Jawab |
|-------|----------------|
| `SequenceRecord` | Menyimpan satu sekuens beserta metadata. Menghitung frekuensi nukleotida (Dictionary) dan GC content. |
| `Parser` | Mengubah berkas FASTA, FASTQ, atau ZIP menjadi List berisi objek SequenceRecord. Menggunakan BioPython bila tersedia, fallback ke parser manual bila tidak. |
| `Analyzer` | Mengurutkan List berdasarkan GC content, mengambil Top-N, menghitung statistik agregat, menghasilkan CSV dan grafik EDA. |

### Struktur Data yang Digunakan

- **Dictionary** untuk menghitung frekuensi nukleotida per sekuens. Akses rata-rata O(1) per basa, total O(n) untuk satu sekuens.
- **List** untuk menyimpan seluruh objek SequenceRecord. Pengurutan menggunakan Timsort bawaan Python dengan kompleksitas O(M log M), M adalah jumlah sekuens.

---

## Cara Menggunakan

### Via Website (disarankan)

Buka https://web-production-4b2df.up.railway.app di browser, unggah berkas FASTA/FASTQ atau klik "Coba dataset demo", lalu unduh hasilnya. Tidak perlu instalasi apapun.

### Via Desktop App

Download `SeqFyre.exe` dari [GitHub Releases](https://github.com/FDZ3006/seqfyre/releases/latest), lalu jalankan langsung dengan double click. Tidak perlu install Python. Mendukung drag-and-drop berkas FASTA/FASTQ dan ekspor hasil ke ZIP.

> **Catatan:** Windows mungkin menampilkan peringatan SmartScreen saat pertama kali dijalankan. Klik **More info** lalu **Run anyway** untuk melanjutkan - ini normal untuk aplikasi yang belum memiliki sertifikat digital.

### Via Lokal

Untuk menjalankan di komputer sendiri:

```bash
git clone https://github.com/FDZ3006/seqfyre.git
cd seqfyre
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python app.py
```

Buka `http://127.0.0.1:5000` di browser.

### Via Command Line

```bash
python cli.py data/dataset_16S_rRNA_SeqFyre.fasta --top 3 --outdir hasil
```

Hasil tersimpan di folder `hasil/`: berkas CSV dan 4 grafik PNG.

---

## Hasil Analisis

Hasil pengurutan sekuens dataset bawaan berdasarkan GC content (menurun):

| Rank | Aksesi | Organisme | GC (%) | Kelompok |
|------|--------|-----------|--------|----------|
| 1 | NR_075056.2 | Aquifex aeolicus | 65.15 | Termofilik |
| 2 | NR_037066.1 | Thermus thermophilus | 63.96 | Termofilik |
| 3 | NR_102775.2 | Thermotoga maritima | 63.54 | Termofilik |
| 4 | NR_115284.2 | Geobacillus stearothermophilus | 59.44 | Termofilik |
| 5 | NR_112116.2 | Bacillus subtilis | 55.10 | Mesofilik |
| 6 | NR_114042.1 | Escherichia coli | 54.92 | Mesofilik |
| 7 | NR_024570.1 | Escherichia coli (U 5/41) | 54.73 | Mesofilik |
| 8 | NR_026078.1 | Pseudomonas aeruginosa | 54.07 | Mesofilik |
| 9 | NR_112081.1 | Herbaspirillum rubrisubalbicans | 53.21 | Mesofilik |
| 10 | NR_118997.2 | Staphylococcus aureus | 51.03 | Mesofilik |

Keempat bakteri termofilik menempati rank 1-4 dan keenam bakteri mesofilik menempati rank 5-10. Pipeline berhasil menstratifikasi organisme sesuai dengan suhu pertumbuhannya.

---

## Eksplorasi Data

### Histogram distribusi GC

![Histogram GC](assets/histogram_gc.png)

Distribusi GC content menunjukkan pola bimodal dengan dua kelompok yang terpisah jelas: bakteri mesofilik di kisaran 51-55% dan bakteri termofilik di kisaran 59-65%.

### GC content per sekuens

![Bar chart GC](assets/barchart_gc.png)

### Komposisi nukleotida

![Komposisi](assets/composition.png)

Bakteri termofilik menunjukkan kandungan G dan C yang lebih tinggi dengan kandungan A yang lebih rendah dibandingkan bakteri mesofilik.

### Panjang sekuens vs GC content

![Scatter](assets/scatter_len_gc.png)

Titik-titik data berkumpul di sekitar 1.450-1.580 bp, mencerminkan homogenitas panjang gen 16S rRNA. Pemisahan pada sumbu GC content sesuai dengan kelompok ekologis masing-masing organisme.

---

## Pengujian

```bash
python -m unittest discover -s tests -v
```

13 pengujian unit mencakup kelas SequenceRecord, Parser, dan Analyzer.

---

---

# English Documentation

# SeqFyre

All core logic is implemented in Python. Reads FASTA or FASTQ files, computes nucleotide frequency, sorts sequences by GC content, and produces visualizations and CSV export. HTML, CSS, and JavaScript in this repository serve only as the web display interface and contain no analysis logic.

**Access directly (no installation required):** https://web-production-4b2df.up.railway.app

**Download desktop app:** [SeqFyre.exe - GitHub Releases](https://github.com/FDZ3006/seqfyre/releases/latest) (Windows, no Python installation needed)

---

**Course:** Bioinformatics Data Structures (BIF1222)  
**Study Program:** Bioinformatics, IPB University  
**Author:** Fajar Dhuha Zafran (G0401241009)

---

## Background

The structural stability of nucleic acids is influenced by base composition. Guanine-Cytosine (G-C) base pairs have three hydrogen bonds and stronger base-stacking forces compared to Adenine-Thymine (A-T) pairs which have only two hydrogen bonds. This difference makes molecules with higher GC content more resistant to elevated temperatures.

Thermophilic bacteria, which live in extreme temperature environments such as hot springs and hydrothermal vents, exploit this mechanism by enriching the G and C content in their structural RNA genes. Conversely, mesophilic bacteria living at moderate temperatures have lower GC content.

The 16S rRNA gene was chosen as the dataset because it is a universal phylogenetic biomarker approximately 1,500 bp in length, conserved enough for cross-species comparison while still reflecting temperature selection pressure on its nucleotide composition.

When this pipeline sorts sequences by GC content in descending order, thermophilic bacteria consistently appear at the top, validating that the pipeline reflects real biological phenomena.

---

## Dataset

File: `data/dataset_16S_rRNA_SeqFyre.fasta`

10 16S rRNA gene sequences from NCBI RefSeq, consisting of 4 thermophilic and 6 mesophilic bacteria.

| Accession | Organism | Group | Optimal Temperature |
|-----------|----------|-------|-------------------|
| NR_075056.2 | Aquifex aeolicus VF5 | Thermophile | ~95°C |
| NR_037066.1 | Thermus thermophilus HB8 | Thermophile | ~70-75°C |
| NR_102775.2 | Thermotoga maritima MSB8 | Thermophile | ~80°C |
| NR_115284.2 | Geobacillus stearothermophilus | Thermophile | ~65°C |
| NR_112116.2 | Bacillus subtilis | Mesophile | ~37°C |
| NR_114042.1 | Escherichia coli | Mesophile | ~37°C |
| NR_024570.1 | Escherichia coli (strain U 5/41) | Mesophile | ~37°C |
| NR_026078.1 | Pseudomonas aeruginosa | Mesophile | ~37°C |
| NR_112081.1 | Herbaspirillum rubrisubalbicans | Mesophile | ~30°C |
| NR_118997.2 | Staphylococcus aureus | Mesophile | ~37°C |

---

## Features
All mandatory points below are fully implemented in Python using List and Dictionary data structures as the core of the pipeline.
Mandatory requirements fulfilled:

- Read FASTA or FASTQ files
- Store data in List structure
- Calculate nucleotide frequency using Dictionary
- Sort sequences by GC content
- Display top 3 sequences (configurable: 3, 5, or 10)
- Visualize results as charts based on GC values
- Write results to CSV file

Additional features:

- Web interface (Flask) accessible directly from browser without installation
- Support for FASTA, FASTQ, and ZIP file uploads
- 4 exploratory data analysis charts (histogram, bar chart, nucleotide composition, scatter plot)
- Download analysis results as ZIP containing CSV and charts
- Object-oriented code with separate classes (SequenceRecord, Parser, Analyzer)
- Bilingual interface (Indonesian and English)
- Three display themes (light, dark, Fyre)

---

## Project Structure

```
seqfyre/
├── seqfyre/               # Core package (OOP)
│   ├── models.py          # SequenceRecord - single sequence model
│   ├── parser.py          # Parser - reads FASTA/FASTQ/ZIP
│   ├── analyzer.py        # Analyzer - sorting, statistics, charts, CSV
│   └── packager.py        # Packages results into ZIP
├── app.py                 # Flask web application
├── cli.py                 # Command-line pipeline
├── templates/index.html   # Web interface
├── static/                # CSS and JavaScript
├── data/                  # 16S rRNA dataset
├── tests/                 # Unit tests
├── assets/                # Charts for documentation
└── requirements.txt
```

### Class Overview

| Class | Responsibility |
|-------|----------------|
| `SequenceRecord` | Stores a single sequence with metadata. Computes nucleotide frequency (Dictionary) and GC content. |
| `Parser` | Converts FASTA, FASTQ, or ZIP files into a List of SequenceRecord objects. Uses BioPython when available, falls back to a manual parser otherwise. |
| `Analyzer` | Sorts the List by GC content, retrieves Top-N sequences, computes aggregate statistics, and produces CSV and EDA charts. |

### Data Structures Used

- **Dictionary** for computing nucleotide frequency per sequence. Average O(1) access per base, total O(n) for one sequence.
- **List** for storing all SequenceRecord objects. Sorted using Python's built-in Timsort with O(M log M) complexity, where M is the number of sequences.

---

## How to Use

### Via Website (recommended)

Open https://web-production-4b2df.up.railway.app in your browser, upload a FASTA/FASTQ file or click "Try demo dataset", then download the results. No installation required.

### Via Desktop App

Download `SeqFyre.exe` from [GitHub Releases](https://github.com/FDZ3006/seqfyre/releases/latest) and run it with a double-click. No Python installation needed. Supports drag-and-drop FASTA/FASTQ files and ZIP export of results.

> **Note:** Windows may show a SmartScreen warning on first run. Click **More info** then **Run anyway** to continue - this is normal for apps without a digital certificate.

### Via Local

To run on your own machine:

```bash
git clone https://github.com/FDZ3006/seqfyre.git
cd seqfyre
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python app.py
```

Open `http://127.0.0.1:5000` in your browser.

### Via Command Line

```bash
python cli.py data/dataset_16S_rRNA_SeqFyre.fasta --top 3 --outdir hasil
```

Results are saved to the `hasil/` folder: a CSV file and 4 PNG charts.

---

## Analysis Results

Sequences from the built-in dataset sorted by GC content (descending):

| Rank | Accession | Organism | GC (%) | Group |
|------|-----------|----------|--------|-------|
| 1 | NR_075056.2 | Aquifex aeolicus | 65.15 | Thermophile |
| 2 | NR_037066.1 | Thermus thermophilus | 63.96 | Thermophile |
| 3 | NR_102775.2 | Thermotoga maritima | 63.54 | Thermophile |
| 4 | NR_115284.2 | Geobacillus stearothermophilus | 59.44 | Thermophile |
| 5 | NR_112116.2 | Bacillus subtilis | 55.10 | Mesophile |
| 6 | NR_114042.1 | Escherichia coli | 54.92 | Mesophile |
| 7 | NR_024570.1 | Escherichia coli (U 5/41) | 54.73 | Mesophile |
| 8 | NR_026078.1 | Pseudomonas aeruginosa | 54.07 | Mesophile |
| 9 | NR_112081.1 | Herbaspirillum rubrisubalbicans | 53.21 | Mesophile |
| 10 | NR_118997.2 | Staphylococcus aureus | 51.03 | Mesophile |

All four thermophilic bacteria occupy ranks 1-4, and all six mesophilic bacteria occupy ranks 5-10. The pipeline successfully stratifies organisms according to their growth temperature.

---

## Exploratory Data Analysis

### GC Distribution Histogram

![Histogram GC](assets/histogram_gc.png)

The GC content distribution shows a bimodal pattern with two clearly separated groups: mesophilic bacteria in the 51-55% range and thermophilic bacteria in the 59-65% range.

### GC Content per Sequence

![Bar chart GC](assets/barchart_gc.png)

### Nucleotide Composition

![Composition](assets/composition.png)

Thermophilic bacteria show higher G and C content with lower A content compared to mesophilic bacteria.

### Sequence Length vs GC Content

![Scatter](assets/scatter_len_gc.png)

Data points cluster around 1,450-1,580 bp, reflecting the length homogeneity of 16S rRNA genes. Separation along the GC content axis corresponds to the ecological group of each organism.

---

## Testing

```bash
python -m unittest discover -s tests -v
```

13 unit tests covering the SequenceRecord, Parser, and Analyzer classes.

---

*SeqFyre - Mini Project Struktur Data Bioinformatika (BIF1222), IPB University*
