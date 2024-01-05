# SacCerML
SacCerML is a ML classifier that enables gene calling for Saccharomyces cerevisiae FASTA sequences

## Gene Prediction in Saccharomyces cerevisiae

### Introduction
This Python script combines bioinformatics tools and machine learning to predict genes in Saccharomyces cerevisiae (baker's yeast). The approach leverages open reading frame (ORF) prediction using the `orfipy_core` library and utilizes a Decision Tree classifier from `scikit-learn` for gene type classification.

### Main objectives
- **Data Processing:** Load genomic sequences from a FASTA file, predict ORFs, and extract various features, including Codon Adaptation Index (CAI), checksum, hydrophobicity, isoelectric point, aromaticity, instability index, molecular weight, secondary structure, and molar extinction coefficients.
- **Model Training and Prediction:** Train a Decision Tree classifier on a labeled dataset (`scerevisiae.csv`) containing features and ORF types. Evaluate the model's accuracy on a test set and optionally save the trained model for future use.

## Methods
1. **Data Loading and Preprocessing:**
   - Read genomic sequences from a FASTA file, predict ORFs, and extract features.
   - Prepare a labeled dataset (`scerevisiae.csv`) for training and testing.

2. **Model Training:**
   - Split the dataset into training and test sets.
   - Utilize a Decision Tree classifier to build and train the gene prediction model.

3. **Prediction and Evaluation:**
   - Predict gene types on the test set and calculate the accuracy of the model.

4. **Optional Model Saving:**
   - Optionally save the trained model using the `joblib` library if the `enable_model_saving` flag is set.

### Results
- Print protein sequences and their ORF type predictions

### Discussion/Conclusion
The script demonstrates an integrated approach to gene prediction, combining ORF analysis with machine learning. It achieves accurate gene type predictions based on a set of genomic features. Considerations for the `longest_M_starting_orf_only` parameter are noted, providing flexibility in the analysis. The option to save the trained model facilitates future use and sharing.



FULL DOCUMENTATION IS BEING BUILT... For now, use [the notebook](./notebook/SacCerML.ipynb) to play around with things;

See it on Google Colab [here](https://colab.research.google.com/drive/107045cvJIsOS30DV3DBoW5PHv6DsdM2F?usp=sharing)
