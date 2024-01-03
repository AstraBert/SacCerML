import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.CheckSum import crc32
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex
from Bio.SeqUtils.CodonUsageIndices import SharpEcoliIndex
from Bio.SeqUtils import six_frame_translations
from Bio.Seq import Seq
from Bio import SeqIO
import gzip
from math import floor
from sklearn.metrics import accuracy_score
from argparse import ArgumentParser

argparse = ArgumentParser()
argparse.add_argument("-s", "--enable_model_saving", help="Save the model to load it in future; keep in mind that the model can take up to 11 GB of disk space!", action="store_true", default=False, required=False)
args = argparse.parse_args()

ems = args.enable_model_saving

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

def process_dna(fasta_file):
    fas = load_data(fasta_file)
    seqs = [seq for seq in list(fas.values())]
    heads = [seq for seq in list(fas.keys())]
    data = []
    for seq in seqs:
        cai = calculate_cai(seq)
        cksm = checksum(seq)
        hydr = hidrophobicity(seq)
        isl = isoelectric_pt(seq)
        arm = aromatic(seq)
        inst = instable(seq)
        mw = weight(seq)
        se_st = sec_struct(seq).split(",")
        se_st1 = se_st[0]
        se_st2 = se_st[1]
        se_st3 = se_st[2]
        me = mol_ext(seq).split(",")
        me1 = me[0]
        me2 = me[1]
        n = pd.DataFrame({"CAI": [cai], "CHECKSUM": [cksm], "HIDROPHOBICITY": [hydr], "ISOELECTRIC": [isl],"AROMATIC": [arm],"INSTABLE": [inst], "MW": [mw], "HELIX": [se_st1], "TURN": [se_st2], "SHEET": [se_st3],"MOL_EXT_RED": [me1], "MOL_EXT_OX": [me2]})
        data.append(n)
    return data, heads

print("Loading data...")
# Load the data from the CSV file
data = pd.read_csv('scerevisiae.csv')
print("Loaded data")

print("Generating training and test data...")
# Features (replace with the actual feature columns in your CSV)
X = data.iloc[:, 1:]

# Labels (replace with the actual label column in your CSV)
y = data['GENE_NAME']

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
print("Generated training and test data")

print("Building and training the model...")
# Create and train the Decision Tree classifier (use the same model you trained before)
model = RandomForestClassifier(n_estimators=25, n_jobs=-1,random_state=42)
model = model.fit(X, y)  # Uncomment this line if clf needs training
print("Built and trained the model... Now predicting")

def predict_genes(infile, model=model):
    X, headers = process_dna(infile)
    predictions = []
    for x in X:
      p = model.predict(x)
      predictions.append(p)
    for i in range(len(predictions)):
        print(f"{headers[i]} is predicted as {predictions[i][0]}")

# Make predictions on the test set
y_pred = model.predict(X_test)

# Evaluate the accuracy of the model
accuracy = accuracy_score(y_test, y_pred)
print(f"Accuracy: {accuracy}")

if ems:
  from joblib import dump
  print("Saving model...")
  dump(model, "SacCerML.joblib")
  print("Saved")

print("All done")
