# conda activate grover

python _grover_in.py

cd ../grover/grover

python scripts/save_features.py --data_path ../../scripts/_grover.csv \
                                --save_path ../../scripts/_grover.npz \
                                --features_generator rdkit_2d_normalized \
                                --restart

python main.py predict --data_path ../../scripts/_grover.csv \
                       --features_path ../../scripts/_grover.npz \
                       --checkpoint_dir data/osm_reg \
                       --no_features_scaling \
                       --output ../../scripts/_grover_reg.csv

python main.py predict --data_path ../../scripts/_grover.csv \
                       --features_path ../../scripts/_grover.npz \
                       --checkpoint_dir data/osm_clf \
                       --no_features_scaling \
                       --output ../../scripts/_grover_clf.csv

python main.py predict --data_path ../../scripts/_grover.csv \
                       --features_path ../../scripts/_grover.npz \
                       --checkpoint_dir data/osm \
                       --no_features_scaling \
                       --output ../../scripts/_grover_ic.csv

cd ../../scripts

python _grover_out.py

rm _grover*.csv
rm _grover*.npz
