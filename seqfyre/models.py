"""
models.py
---------
Definisi kelas SequenceRecord.
Satu objek merepresentasikan satu sekuens nukleotida beserta metadata
dan perhitungan turunannya (frekuensi nukleotida, GC content).
"""


class SequenceRecord:
    # Basa nukleotida standar yang dihitung.
    # N dipisah agar tidak merusak akurasi GC content.
    BASES = ("A", "T", "G", "C")

    def __init__(self, seq_id, description, sequence):
        self.id = seq_id
        self.description = description
        # Normalisasi: huruf besar, buang whitespace
        self.sequence = "".join(sequence.split()).upper()
        self.organism = self._extract_organism(description)

    def _extract_organism(self, description):
        # Ambil 2 kata setelah accession sebagai nama organisme
        # Contoh header NCBI: 'NR_037066.1 Thermus thermophilus HB8 16S ...'
        parts = description.split()
        if len(parts) >= 3:
            return parts[1] + " " + parts[2]
        if len(parts) == 2:
            return parts[1]
        return "Unknown organism"

    def get_length(self):
        # Panjang total sekuens termasuk basa ambigu N
        return len(self.sequence)

    def nucleotide_frequency(self):
        # Hitung frekuensi tiap nukleotida menggunakan Dictionary
        # Kompleksitas O(n), akses Dictionary rata-rata O(1) per basa
        freq = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0}
        for base in self.sequence:
            if base in freq:
                freq[base] += 1
            else:
                # Basa IUPAC lain dianggap ambigu
                freq["N"] += 1
        return freq

    def gc_content(self):
        # Persentase GC = (G + C) / (A + T + G + C) * 100
        # Basa N dikecualikan dari penyebut
        freq = self.nucleotide_frequency()
        valid = freq["A"] + freq["T"] + freq["G"] + freq["C"]
        if valid == 0:
            return 0.0
        return (freq["G"] + freq["C"]) / valid * 100

    def at_content(self):
        # Persentase AT, komplemen dari GC content
        gc = self.gc_content()
        return 100.0 - gc if gc else 0.0

    def composition_percent(self):
        # Persentase tiap basa terhadap panjang sekuens
        length = self.get_length()
        if length == 0:
            return {"A": 0.0, "T": 0.0, "G": 0.0, "C": 0.0, "N": 0.0}
        freq = self.nucleotide_frequency()
        result = {}
        for b in freq:
            result[b] = freq[b] / length * 100
        return result

    def classify_by_gc(self, threshold=57.0):
        # Klasifikasi ekologis berdasarkan ambang GC
        # Termofil cenderung GC tinggi (>57%), mesofil lebih rendah
        # Ambang 57% dipilih karena berada di celah bimodal dataset
        if self.gc_content() >= threshold:
            return "Termofil"
        return "Mesofil"

    def to_row(self, rank=None):
        # Ubah objek menjadi satu baris dictionary untuk tabel dan CSV
        freq = self.nucleotide_frequency()
        row = {
            "Rank": rank if rank is not None else "",
            "Accession_ID": self.id,
            "Organism": self.organism,
            "Length": self.get_length(),
            "GC_Content(%)": round(self.gc_content(), 2),
            "A_Count": freq["A"],
            "T_Count": freq["T"],
            "G_Count": freq["G"],
            "C_Count": freq["C"],
            "N_Count": freq["N"],
            "Classification": self.classify_by_gc(),
        }
        return row

    def __len__(self):
        return self.get_length()

    def __repr__(self):
        return (
            "SequenceRecord(id=" + repr(self.id) +
            ", organism=" + repr(self.organism) +
            ", len=" + str(self.get_length()) +
            ", gc=" + str(round(self.gc_content(), 2)) + "%)"
        )
