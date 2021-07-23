import csv
import pandas as pd

score_clf = []
with open("_grover_clf.csv", "r") as f:
    reader = csv.reader(f)
    next(reader)
    for r in reader:
        score_clf += [float(r[1])]

score_reg = []
with open("_grover_reg.csv", "r") as f:
    reader = csv.reader(f)
    next(reader)
    for r in reader:
        score_reg += [float(r[1])]

score_ic = []
with open("_grover_ic.csv", "r") as f:
    reader = csv.reader(f)
    next(reader)
    for r in reader:
        score_ic += [float(r[1])]


df_ = pd.DataFrame({
    "ActivityClfVocab": score_clf,
    "ActivityRegVocab": score_reg,
    "ActivityIcVocab": score_ic
    }
)

df_.to_csv("_pred_grover.csv", index=False)
