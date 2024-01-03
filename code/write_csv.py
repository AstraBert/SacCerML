from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.CheckSum import crc32
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex
from Bio.SeqUtils.CodonUsageIndices import SharpEcoliIndex
from Bio.SeqUtils import six_frame_translations
from Bio.Seq import Seq
import gzip
from Bio import SeqIO
from math import floor


## CSV structure: GENE_NAME,CAI,CHECKSUM,HIDROPHOBICITY,ISOELECTRIC,AROMATIC,INSTABLE,MW,HELIX,TURN,SHEET,MOL_EXT_RED,MOL_EXT_OX
def load_data(infile):
    """Load data from infile if it is in fasta format (after having unzipped it, if it is zipped)"""
    if infile.endswith(".gz"):  # If file is gzipped, unzip it
        y = gzip.open(infile, "rt", encoding="latin-1")
        # Read file as fasta if it is fasta
        if infile.endswith(".fasta.gz") or infile.endswith(".fna.gz") or infile.endswith(".fas.gz") or infile.endswith(".fa.gz"):
            records = SeqIO.parse(y, "fasta")
            sequences = {}
            for record in records:
                sequences.update({str(record.id): str(record.seq)})
            y.close()
            return sequences
        else:
            y.close()
            raise ValueError("File is the wrong format")
    # Read file directly as fasta if it is a not zipped fasta: handle also more uncommon extensions :-)
    elif infile.endswith(".fasta") or infile.endswith(".fna") or infile.endswith(".fas") or infile.endswith(".fa"):
        with open(infile, "r") as y:
            records = SeqIO.parse(y, "fasta")
            sequences = {}
            for record in records:
                sequences.update({str(record.id): str(record.seq)})
            y.close()
            return sequences
    else:
        raise ValueError("File is the wrong format")



def calculate_cai(dna,index=SharpEcoliIndex):
    cai = CodonAdaptationIndex()
    cai.set_cai_index(index)
    if len(dna) % 3 == 0:
        a = cai.cai_for_gene(dna)
    else:
        six_translated = six_frame_translations(dna)
        n = six_translated.split("\n")
        frames = {"0;F": n[5], "1;F": n[6], "2;F": n[7],"0;R": n[12], "1;R": n[11], "2;R": n[10]}
        ind = 0
        for i in list(frames.keys()):
            k = frames[i].replace(" ","")
            if "M" in k and "*" in k:
                if i.split(";")[0]=="F" and k.index("M") < k.index("*"):
                    if len(k) <= len(dna)/3:
                        ind = int(i.split("")[0])
                        break
                elif i.split(";")[0]=="R" and k.index("M") > k.index("*"):
                    if len(k) <= len(dna)/3:
                        ind = len(dna) - int(i.split("")[0])
                        break
        if ind == 0:
            cods = 3*floor(len(dna)/3)
            dna = dna[:cods]
            a = cai.cai_for_gene(dna)
        elif 1 <= ind <= 2:
            if len(dna[ind:])%3==0:
                dna = dna[ind:]
            else:
                cods = 3*floor((len(dna)-ind)/3)
                dna = dna[ind:cods+ind]
                a = cai.cai_for_gene(dna)
        else:
            if len(dna[:ind])%3==0:
                dna = dna[ind:]
            else:
                cods = 3*floor((len(dna)-ind)/3)
                dna = dna[:cods]
                a = cai.cai_for_gene(dna)
    return a


def checksum(dna):
    return crc32(dna)


def hidrophobicity(dna):
    protein_sequence = str(Seq(dna).translate())
    protein_sequence = protein_sequence.replace("*","")
    hydrophobicity_score = ProteinAnalysis(protein_sequence).gravy()
    return hydrophobicity_score

def isoelectric_pt(dna):
    protein_sequence = str(Seq(dna).translate())
    protein_sequence = protein_sequence.replace("*","")
    isoelectric = ProteinAnalysis(protein_sequence).isoelectric_point()
    return isoelectric

def aromatic(dna):
    protein_sequence = str(Seq(dna).translate())
    protein_sequence = protein_sequence.replace("*","")
    arom = ProteinAnalysis(protein_sequence).aromaticity()
    return arom


def instable(dna):
    protein_sequence = str(Seq(dna).translate())
    protein_sequence = protein_sequence.replace("*","")
    inst = ProteinAnalysis(protein_sequence).instability_index()
    return inst

def weight(dna):
    protein_sequence = str(Seq(dna).translate())
    protein_sequence = protein_sequence.replace("*","")
    wgt = ProteinAnalysis(protein_sequence).molecular_weight()
    return wgt

def sec_struct(dna):
    protein_sequence = str(Seq(dna).translate())
    protein_sequence = protein_sequence.replace("*","")
    second_struct = ProteinAnalysis(protein_sequence).secondary_structure_fraction()
    return ",".join([str(s) for s in second_struct])

def mol_ext(dna):
    protein_sequence = str(Seq(dna).translate())
    protein_sequence = protein_sequence.replace("*","")
    molar_ext = ProteinAnalysis(protein_sequence).molar_extinction_coefficient()
    return ",".join([str(s) for s in molar_ext])


if __name__=="__main__":
    print("Loading data...")
    fasta = load_data("data/s_cerevisiae.fasta")
    print(list(fasta.keys())[0])
    print("Done")
    print("Writing csv...")
    with open("data/scerevisiae.csv","w") as csv:
        csv.write("ORF_TYPE,CAI,CHECKSUM,HIDROPHOBICITY,ISOELECTRIC,AROMATIC,INSTABLE,MW,HELIX,TURN,SHEET,MOL_EXT_RED,MOL_EXT_OX\n")
        c = 0
        for i in list(fasta.keys()):
            c+=1
            orf = i
            cai = calculate_cai(fasta[i])
            cksm = checksum(fasta[i])
            hydr = hidrophobicity(fasta[i])
            isl = isoelectric_pt(fasta[i])
            arm = aromatic(fasta[i])
            inst = instable(fasta[i])
            mw = weight(fasta[i])
            se_st = sec_struct(fasta[i])
            me = mol_ext(fasta[i])
            csv.write(f"{orf},{cai},{cksm},{hydr},{isl},{arm},{inst},{mw},{se_st},{me}\n")
            if c%500==0:
                print(f"Processed {c} reads")
    csv.close()
    print("Done")



