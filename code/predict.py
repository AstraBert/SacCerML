from joblib import load
from model import process_dna, model
from argparse import ArgumentParser

argparse = ArgumentParser()
argparse.add_argument(
    "-i",
    "--infile",
    help="Path to input file with S. cerevisiae sequences",
    required=True,
)
args = argparse.parse_args()

inf = args.infile


# loaded_model = load("../SacCerML.joblib")
def predict_genes(infile, model=model):
    X, proteins = process_dna(infile)
    headers = list(X.keys())
    predictions = []
    for x in list(X.values()):
        p = model.predict(x)
        predictions.append(p)
    for i in range(len(predictions)):
        print(
            f"{headers[i]} protein sequence is\n{proteins[headers[i]]}\nand is predicted as {predictions[i][0]}\n"
        )


if __name__ == "__main__":
    predict_genes(inf)
