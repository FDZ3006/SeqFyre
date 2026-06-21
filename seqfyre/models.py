"""
models.py
=========
Definisi kelas SequenceRecord: representasi satu rekaman sekuens biologis.

Kelas ini bertanggung jawab penuh atas SATU sekuens — menyimpan metadata
(accession, organisme, deskripsi) dan menghitung properti turunannya:
frekuensi nukleotida (struktur Dictionary), GC content, dan komposisi basa.

Desain OOP: satu objek = satu sekuens. Analyzer-lah yang nanti mengelola
kumpulan objek ini dalam sebuah List.
"""

from __future__ import annotations


class SequenceRecord:
    """Merepresentasikan satu sekuens 16S rRNA (atau sekuens nukleotida apa pun).

    Attributes:
        id (str): Nomor aksesi / SeqID dari header FASTA (mis. 'NR_037066.1').
        description (str): Header lengkap setelah tanda '>'.
        organism (str): Nama organisme hasil ekstraksi dari header.
        sequence (str): Untaian nukleotida (huruf besar, tanpa whitespace).
    """

    # Basa nukleotida standar yang dihitung. 'N' (ambigu, IUPAC) sengaja
    # dipisah agar bisa dilaporkan namun tidak merusak akurasi GC content.
    BASES = ("A", "T", "G", "C")

    def __init__(self, seq_id: str, description: str, sequence: str) -> None:
        self.id: str = seq_id
        self.description: str = description
        # Normalisasi: huruf besar + buang karakter non-huruf (newline, spasi).
        self.sequence: str = "".join(sequence.split()).upper()
        self.organism: str = self._extract_organism(description)

    # ------------------------------------------------------------------ #
    # Helper privat
    # ------------------------------------------------------------------ #
    @staticmethod
    def _extract_organism(description: str) -> str:
        """Ambil nama organisme dari header FASTA.

        Header NCBI berbentuk:
            'NR_037066.1 Thermus thermophilus HB8 16S ribosomal RNA, ...'
        Kita ambil 2 kata setelah accession sebagai 'Genus species'.
        """
        parts = description.split()
        if len(parts) >= 3:
            return f"{parts[1]} {parts[2]}"
        if len(parts) == 2:
            return parts[1]
        return "Unknown organism"

    # ------------------------------------------------------------------ #
    # Properti turunan
    # ------------------------------------------------------------------ #
    @property
    def length(self) -> int:
        """Panjang total sekuens (termasuk basa ambigu 'N')."""
        return len(self.sequence)

    def nucleotide_frequency(self) -> dict[str, int]:
        """Hitung frekuensi tiap nukleotida memakai struktur Dictionary.

        Inilah implementasi inti persyaratan WAJIB: 'menghitung frekuensi
        nukleotida menggunakan Dictionary'. Kompleksitas waktu O(n) dengan
        akses hash table rata-rata O(1) per basa.

        Returns:
            dict: mis. {'A': 310, 'T': 249, 'G': 525, 'C': 462, 'N': 0}
        """
        freq: dict[str, int] = {base: 0 for base in self.BASES}
        freq["N"] = 0  # penampung basa ambigu / non-standar
        for base in self.sequence:
            if base in freq:
                freq[base] += 1
            else:
                # Basa IUPAC lain (R, Y, K, M, ...) diperlakukan sebagai ambigu.
                freq["N"] += 1
        return freq

    @property
    def gc_content(self) -> float:
        """Persentase GC = (G + C) / (A + T + G + C) * 100.

        Basa 'N' DIKECUALIKAN dari penyebut agar persentase GC murni
        mencerminkan termodinamika pasangan basa, bukan noise sekuensing.
        """
        freq = self.nucleotide_frequency()
        valid = freq["A"] + freq["T"] + freq["G"] + freq["C"]
        if valid == 0:
            return 0.0
        return (freq["G"] + freq["C"]) / valid * 100

    @property
    def at_content(self) -> float:
        """Persentase AT (komplemen dari GC content)."""
        gc = self.gc_content
        return 100.0 - gc if gc else 0.0

    def composition_percent(self) -> dict[str, float]:
        """Persentase tiap basa relatif terhadap panjang sekuens."""
        if self.length == 0:
            return {b: 0.0 for b in (*self.BASES, "N")}
        freq = self.nucleotide_frequency()
        return {b: freq[b] / self.length * 100 for b in freq}

    def classify_by_gc(self, threshold: float = 57.0) -> str:
        """Klasifikasi ekologis sederhana berdasarkan ambang GC.

        Termofil cenderung GC tinggi (>57%), mesofil lebih rendah.
        Ambang 57% dipilih karena berada di celah bimodal dataset.
        """
        return "Termofil" if self.gc_content >= threshold else "Mesofil"

    def to_row(self, rank: int | None = None) -> dict:
        """Ubah objek menjadi satu baris (dict) untuk tabel / CSV."""
        freq = self.nucleotide_frequency()
        row = {
            "Rank": rank if rank is not None else "",
            "Accession_ID": self.id,
            "Organism": self.organism,
            "Length": self.length,
            "GC_Content(%)": round(self.gc_content, 2),
            "A_Count": freq["A"],
            "T_Count": freq["T"],
            "G_Count": freq["G"],
            "C_Count": freq["C"],
            "N_Count": freq["N"],
            "Classification": self.classify_by_gc(),
        }
        return row

    # ------------------------------------------------------------------ #
    # Dunder methods (membuat objek nyaman dipakai)
    # ------------------------------------------------------------------ #
    def __len__(self) -> int:
        return self.length

    def __repr__(self) -> str:
        return (
            f"SequenceRecord(id={self.id!r}, organism={self.organism!r}, "
            f"len={self.length}, gc={self.gc_content:.2f}%)"
        )
