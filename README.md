<div align="center">

# 🔥 SeqFyre

**Pipeline analisis struktur data untuk gen 16S rRNA — mengurutkan sekuens berdasarkan GC content dan menyingkap gradien termal bakteri.**

Mini Project · Struktur Data Bioinformatika (BIF1223) · IPB University

[Fitur](#-fitur) · [Demo](#-cara-menjalankan) · [Arsitektur](#-arsitektur--struktur-proyek) · [Deploy](#-deploy-ke-render) · [EDA](#-eksplorasi-data-eda)

</div>

---

## 📖 Ringkasan

**SeqFyre** (dari *Sequence* + *Fyre*/api) adalah aplikasi analisis bioinformatika yang membaca berkas FASTA/FASTQ gen 16S rRNA, menghitung frekuensi nukleotida, mengurutkan sekuens berdasarkan **GC content**, dan memvisualisasikan hasilnya. Tersedia dalam dua antarmuka:

- **CLI** (`cli.py`) — pipeline mandiri yang memenuhi seluruh persyaratan dasar tugas.
- **Web app** (`app.py`) — antarmuka Flask yang bisa di-deploy ke Render; pengguna cukup membuka tautan, mengunggah berkas, dan mengunduh hasil tanpa instalasi apa pun.

Inti dari proyek ini bukan sekadar manipulasi string: ketika sekuens diurutkan menurun berdasarkan GC content, **bakteri termofilik secara deterministik naik ke puncak** sementara mesofilik mengendap di dasar — pipeline komputasi memvalidasi fenomena adaptasi termal yang nyata di alam.

---

## 🧬 Data Understanding — mengapa dataset ini

Stabilitas struktural asam nukleat sangat dipengaruhi suhu lingkungan. Pasangan basa **Guanin–Sitosin (G–C)** terikat oleh tiga ikatan hidrogen dan memiliki gaya *base-stacking* yang kuat, sehingga jauh lebih stabil secara termodinamika dibanding pasangan **Adenin–Timin (A–T)** yang hanya memiliki dua ikatan hidrogen. Akibatnya, bakteri **termofilik** — yang hidup di sumber air panas dan ventilasi hidrotermal — berevolusi memperkaya kandungan G dan C pada gen struktural RNA-nya untuk mencegah denaturasi pada suhu ekstrem.

Gen **16S rRNA** dipilih karena merupakan biomarker filogenetik universal (~1.500 bp), cukup terkonservasi namun tetap memuat variasi yang mencerminkan tekanan seleksi suhu. Dataset menggabungkan empat termofil dan enam mesofil sehingga pengurutan GC menghasilkan **distribusi bimodal** yang khas — dua puncak terpisah, bukan satu kurva normal.

> Konteks biologis di atas membuat setiap angka GC menjadi bermakna: ia adalah rekam jejak adaptasi suhu yang terekam secara molekuler, bukan statistik acak.

### Dataset

Berkas: [`data/dataset_16S_rRNA_SeqFyre.fasta`](data/dataset_16S_rRNA_SeqFyre.fasta) — 10 sekuens RefSeq (`NR_`) gen 16S rRNA, diverifikasi manual di NCBI Nucleotide.

| Kelompok | Organisme | Aksesi | Suhu optimal |
|----------|-----------|--------|--------------|
| Termofil | *Aquifex aeolicus* VF5 | NR_075056.2 | ~95 °C (hipertermofil) |
| Termofil | *Thermus thermophilus* HB8 | NR_037066.1 | ~70–75 °C |
| Termofil | *Thermotoga maritima* MSB8 | NR_102775.2 | ~80 °C (hipertermofil) |
| Termofil | *Geobacillus stearothermophilus* | NR_115284.2 | ~65 °C |
| Mesofil | *Bacillus subtilis* | NR_112116.2 | ~37 °C |
| Mesofil | *Escherichia coli* | NR_114042.1 | ~37 °C |
| Mesofil | *Escherichia coli* (strain U 5/41) | NR_024570.1 | ~37 °C |
| Mesofil | *Pseudomonas aeruginosa* | NR_026078.1 | ~37 °C |
| Mesofil | *Herbaspirillum rubrisubalbicans* | NR_112081.1 | ~30 °C |
| Mesofil | *Staphylococcus aureus* | NR_118997.2 | ~37 °C |

> **Cara mengunduh dataset sendiri dari NCBI:** buka NCBI Nucleotide, masukkan daftar aksesi dipisahkan koma, lalu **Send to → File → Format: FASTA → Create File**. (Catatan: basis data PopSet sudah dipensiunkan sejak Januari 2025; unduhan kini lewat Nucleotide biasa.)

---

## ✨ Fitur

**Wajib (terpenuhi penuh):**
- ✅ Membaca berkas **FASTA** atau **FASTQ**
- ✅ Menyimpan data dalam struktur **List**
- ✅ Menghitung frekuensi nukleotida memakai **Dictionary**
- ✅ Mengurutkan sekuens berdasarkan **GC content**
- ✅ Menampilkan **3 sekuens terbaik**
- ✅ **Visualisasi grafik** berdasarkan nilai GC
- ✅ Menulis hasil ke berkas **CSV**

**Bonus:**
- 🌐 **Web app Flask** siap deploy ke Render (gratis, tanpa instalasi bagi pengguna)
- 📦 Dukungan unggah **FASTA, FASTQ, dan ZIP** (ZIP berisi banyak berkas otomatis digabung)
- 🔢 Tabel **Top-N konfigurabel** (3 / 5 / 10)
- 📊 **4 grafik EDA** (histogram, bar chart, komposisi, scatter)
- 💾 **Unduh ZIP** hasil (CSV + grafik) — dibuat *in-memory*
- 🧱 **Kode OOP** dengan kelas terpisah (`SequenceRecord`, `Parser`, `Analyzer`)
- 🌏 **Dwibahasa** — Indonesia utama, toggle Inggris di web app
- 🧪 **BioPython** sebagai parser utama (dengan fallback pure-Python)
- ✔️ **Suite pengujian** unit (`tests/`)

---

## 🏗 Arsitektur & struktur proyek

```
seqfyre/
├── seqfyre/                  # Paket inti (OOP)
│   ├── __init__.py
│   ├── models.py             # SequenceRecord — satu sekuens (freq, GC content)
│   ├── parser.py             # Parser — baca FASTA/FASTQ/ZIP (BioPython + fallback)
│   ├── analyzer.py           # Analyzer — sort GC, statistik, CSV, 4 grafik EDA
│   └── packager.py           # build_result_zip — kemas hasil jadi ZIP in-memory
├── app.py                    # Web app Flask
├── cli.py                    # Pipeline CLI (memenuhi semua poin wajib)
├── templates/index.html      # Antarmuka web
├── static/
│   ├── style.css             # Tema "gradien termal"
│   └── app.js                # Logika frontend + toggle bahasa
├── data/                     # Dataset 16S rRNA
├── tests/test_pipeline.py    # 13 pengujian unit
├── assets/                   # Gambar EDA untuk dokumentasi
├── requirements.txt
├── render.yaml · Procfile · runtime.txt   # konfigurasi deploy
└── README.md
```

### Tanggung jawab tiap kelas

| Kelas | Peran |
|-------|-------|
| `SequenceRecord` | Menyimpan satu sekuens + metadata; menghitung frekuensi nukleotida (**Dictionary**), GC content, komposisi, dan klasifikasi termal. |
| `Parser` | Mengubah berkas/teks FASTA·FASTQ·ZIP menjadi **List** `SequenceRecord`. Memakai BioPython bila tersedia, jika tidak memakai parser manual. |
| `Analyzer` | Mengurutkan **List** berdasarkan GC, mengambil Top-N, menghitung statistik agregat, menghasilkan CSV dan 4 grafik EDA (base64/PNG). |

### Struktur data & kompleksitas

- **Dictionary** (`{'A':.., 'T':.., 'G':.., 'C':.., 'N':..}`) untuk frekuensi nukleotida — akses/penambahan rata-rata **O(1)** per basa, total **O(n)** untuk satu sekuens.
- **List** untuk menampung seluruh objek `SequenceRecord` — pengurutan memakai Timsort bawaan Python, **O(M log M)** dengan *M* = jumlah sekuens.
- Top-N diperoleh lewat *array slicing* `sorted[:n]`.

---

## 🚀 Cara menjalankan

### 1. Lokal — siapkan lingkungan

```bash
git clone <URL-REPO-ANDA>
cd seqfyre
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

### 2. Jalankan CLI (memenuhi semua poin wajib)

```bash
python cli.py data/dataset_16S_rRNA_SeqFyre.fasta --top 3 --outdir hasil
```

Keluaran: `hasil/hasil_analisis.csv` dan `hasil/grafik/*.png`.

### 3. Jalankan web app secara lokal

```bash
python app.py
# buka http://127.0.0.1:5000
```

Unggah berkas FASTA/FASTQ/ZIP atau klik **Coba dataset demo**, atur Top-N, lalu unduh ZIP hasil.

### 4. Jalankan pengujian

```bash
python -m unittest discover -s tests -v
```

---

## ☁️ Deploy ke Render

SeqFyre dirancang untuk Render free tier. Dua cara:

**A. Blueprint (otomatis, disarankan)**
1. Push repo ini ke GitHub.
2. Di Render: **New → Blueprint**, pilih repo. Render membaca [`render.yaml`](render.yaml) dan menyiapkan semuanya.

**B. Manual**
1. **New → Web Service**, hubungkan repo.
2. **Build Command:** `pip install -r requirements.txt`
3. **Start Command:** `gunicorn app:app --bind 0.0.0.0:$PORT --workers 2 --timeout 120`
4. **Health Check Path:** `/health`

> Semua keluaran (CSV, grafik, ZIP) dibuat **in-memory**, sesuai sifat *ephemeral filesystem* Render. Grafik dikirim ke frontend sebagai PNG base64.

Setelah live, perbarui tautan demo di sini:

```
🔗 Demo: https://<nama-app>.onrender.com
```

---

## 📊 Eksplorasi Data (EDA)

Keempat grafik berikut dihasilkan otomatis dari dataset bawaan.

### 1. Histogram distribusi GC — pola bimodal
Dua puncak terpisah: mesofil mengelompok di ~51–55%, termofil di ~59–65%. Inilah tanda tangan dataset berkonteks.

![Histogram GC](assets/histogram_gc.png)

### 2. GC content per sekuens
Bar oranye (termofil) merentang melewati ambang yang tak tercapai bar biru (mesofil).

![Bar chart GC](assets/barchart_gc.png)

### 3. Komposisi nukleotida
Pergeseran komposisi terlihat jelas — termofil menekan A dan memperkaya G+C.

![Komposisi nukleotida](assets/composition.png)

### 4. Panjang vs GC content
Titik-titik berkumpul vertikal di ~1.500 bp (homogenitas panjang 16S rRNA), tersebar pada sumbu GC sesuai kelas termal — sekaligus alat deteksi outlier/QC.

![Scatter panjang vs GC](assets/scatter_len_gc.png)

---

## 🏆 Hasil

Pengurutan menurun berdasarkan GC content menempatkan tiga hipertermofil/termofil di puncak:

| Rank | Aksesi | Organisme | GC (%) | Kelas |
|------|--------|-----------|--------|-------|
| 1 | NR_075056.2 | *Aquifex aeolicus* | 65.15 | Termofil |
| 2 | NR_037066.1 | *Thermus thermophilus* | 63.96 | Termofil |
| 3 | NR_102775.2 | *Thermotoga maritima* | 63.54 | Termofil |

Keempat termofil (termasuk *Geobacillus*, 59.44%) menempati rank 1–4, dan keenam mesofil rank 5–10 — pipeline berhasil menstratifikasi organisme berdasarkan suhu pertumbuhannya.

---

## 🛠 Teknologi

Python 3.12 · Flask · BioPython · Matplotlib · Gunicorn · HTML/CSS/JavaScript (vanilla)

## 📄 Lisensi

[MIT](LICENSE) © 2026 SeqFyre — dikembangkan untuk Mini Project BIF1223, IPB University.
