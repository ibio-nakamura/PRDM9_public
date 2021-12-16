import numpy as np
from sklearn.model_selection import train_test_split
import gzip as gz
from core import MCRCDropout, CustomSumPool
from collections import OrderedDict

inverted = {"A":"T", "T":"A", "C":"G", "G":"C"}
padding = [0.0,0,0,0]
ix = ["A", "C", "G","T"]
one_hot_conv = {"A": [1,0,0,0],
           "T": [0,0,0,1],
           "C": [0,1,0,0],
           "G": [0,0,1,0]}

def reverse_complement(sequence):
    return [inverted[base] for base in reversed(sequence)]

def one_hot(sequence):
    return np.array([np.array(one_hot_conv[base], dtype=np.float) for base in sequence] )
    
from keras.models import model_from_json, model_from_yaml
import yaml

def save_model(model, model_name):
    with open(model_name+".json", 'w') as j_file:
        j_file.write(model.to_json())

    model.save_weights(model_name+".h5")
    
def save_model_yaml(model, model_name):
    with open(model_name+".yaml", 'w') as j_file:
        j_file.write(model.to_yaml())

    model.save_weights(model_name+".h5")
    
def load_model(model_name):
    
    with open(model_name + ".json", 'r') as j_file:
        loaded_model_json = j_file.read()
    loaded_model = model_from_json(loaded_model_json)
    loaded_model.load_weights(model_name + ".h5")
    
    return loaded_model

#エラー出たため、cust_objectsのところは{}だったのを自分でMCRCDropoutとCustomSumPool書いた。
#参考：https://blog.shikoan.com/keras-load-custom-loss-model/
def load_model_yaml(model_name, cust_objects={"MCRCDropout":MCRCDropout, "CustomSumPool":CustomSumPool}):
    with open(model_name + ".yaml", 'r') as y_file:
        loaded_model_yaml = y_file.read()
    
    loaded_model = model_from_yaml(loaded_model_yaml,custom_objects=cust_objects)
    loaded_model.load_weights(model_name + ".h5")
    
    return loaded_model

def train_test_val_split(x, y, rs=1):

    X_train, X_val, y_train, y_val = train_test_split(x, y, test_size=0.15, random_state=rs)

    X_train, X_test, y_train, y_test = train_test_split(X_train, y_train, test_size=0.1, random_state=rs)
    
    return X_train, X_val, X_test, y_train, y_val, y_test


def load_fasta_gz(f_name,input_len):
    sequences = []
    cur_string = ""
    s = 0
    with gz.open(f_name) as fasta_file:
        for line in fasta_file:
            line = line.decode("ascii")
            if line[0] == '>':
                s+=1
                if cur_string:
                    assert len(cur_string) ==input_len
                    sequences.append(cur_string)

                cur_string = ""
            else:
                line = line.strip()
                cur_string += line

        assert len(cur_string) ==input_len
        sequences.append(cur_string)


    return sequences


def load_fasta_gz_and_return_OrderedDict(f_name,input_len):
    #sequences = []
    ret_od = OrderedDict()
    key_list = []
    cur_string = ""
    s = 0
    with gz.open(f_name) as fasta_file:
        for line in fasta_file:
            line = line.decode("ascii")
            if line[0] == '>':
                key = line.strip().lstrip('>')
                key_list.append(key)
                if cur_string:
                    assert len(cur_string) ==input_len
                    #sequences.append(cur_string)
                    current_key = key_list[s]
                    ret_od[current_key] = cur_string
                    s+=1
                cur_string = ""
            else:
                line = line.strip()
                cur_string += line
        
        #一番最後の処理
        current_key = key_list[s]
        assert len(cur_string) ==input_len
        ret_od[current_key] = cur_string
        #sequences.append(cur_string)


    return ret_od


def load_recomb_data(aug, input_len):
    hot_sequences = []
    cold_sequences = []
    #for chr_index in [str(1),str(2),str(3),str(4),str(5),str(6),str(7),str(8),str(9),str(10),str(11),str(12),str(13),str(14),str(15),str(16),str(17),str(18),str(19),str(20),str(21),str(22),'X']:
    for chr_index in [str(22)]:
        current_chr_hot_sequences = load_fasta_gz(f_name="../100bp_CNN/100bp_seqs/100bp_train_test_fastas_dir/100bp_hotspots_train_original_chr"+chr_index+".fasta.gz", \
                                                  input_len=input_len)
        hot_sequences += current_chr_hot_sequences

        current_chr_cold_sequences = load_fasta_gz(f_name="../100bp_CNN/100bp_seqs/100bp_train_test_fastas_dir/100bp_coldspots_train_original_chr"+chr_index+".fasta.gz", \
                                                   input_len=input_len)
        cold_sequences += current_chr_cold_sequences

    print("N of hot seqs: ", len(hot_sequences))
    print("N of cold seqs: ", len(cold_sequences))

    if aug:
        X = np.array([one_hot(seq) for seq in hot_sequences] + [one_hot(reverse_complement(seq)) for seq in hot_sequences] + [one_hot(seq) for seq in cold_sequences] + [one_hot(reverse_complement(seq)) for seq in cold_sequences])[:,:input_len] # to make the pooling symmetric
        Y = np.array([1]*(2*len(hot_sequences)) + [0]*(2*len(cold_sequences)))
    else:
        X = np.array([one_hot(seq) for seq in hot_sequences] + [one_hot(seq) for seq in cold_sequences])[:,:input_len] # to make the pooling symmetric
        Y = np.array([1]*len(hot_sequences) + [0]*len(cold_sequences))


    X_train, X_val, X_test, Y_train, Y_val, Y_test = train_test_val_split(X, Y)

    return X_train, X_val, X_test, Y_train, Y_val, Y_test


def predict_mc(self, X_pred, n_preds=100):
    return np.mean([self.predict(X_pred) for i in range(n_preds)], axis=0)
