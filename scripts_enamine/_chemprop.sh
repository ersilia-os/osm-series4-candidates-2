# conda activate chemprop

python _chemprop_in.py

chemprop_predict --test_path _chemprop.csv \
                 --checkpoint_dir ../chemprop/clf_checkpoints \
                 --preds_path _chemprop_clf.csv

chemprop_predict --test_path _chemprop.csv \
                 --checkpoint_dir ../chemprop/reg_checkpoints \
                 --preds_path _chemprop_reg.csv

python _chemprop_out.py

rm _chemprop*.csv
