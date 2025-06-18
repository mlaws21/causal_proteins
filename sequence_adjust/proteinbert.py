
# this is hard code path from main
from sequence_adjust.protein_bert_helper.existing_model_loading import load_pretrained_model
from sequence_adjust.protein_bert_helper.conv_and_global_attention_model import get_model_with_hidden_layers_as_outputs
from scipy.spatial.distance import cosine

# note it can be shorter -- so nonsense is ok?
def protein_bert_scores(ref, seqs, seq_len=253):

    pretrained_model_generator, input_encoder = load_pretrained_model()
    
    seqs = [ref, *seqs]
    
    seq_len += 2
    model = get_model_with_hidden_layers_as_outputs(pretrained_model_generator.create_model(seq_len))
    encoded_x = input_encoder.encode_X(seqs, seq_len)
    _, global_representations = model.predict(encoded_x, batch_size=32, verbose=0)
    dists = [cosine(global_representations[0], global_representations[x]) for x in range(1, len(global_representations))]

    return dists
    
def main():
    ref = "MANLGCWMLVLFVATWSDLGLCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHSQWNKPSKPKTNMKHMAGAAAAGAVVGGLGGYMLGSAMSSPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNITIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRGSSMVLFSSPPVILLISFLIFLIVG"
    seqs = ["MANLGCWMLVLFVATWSDLGLCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHSQWNKPSKPKTNMKHMAGAAAAGAVVGGLGGYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNITIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRGSSMVLFSSPPVILLISFLIFLIVG", "MANLGCWMLVLFVATWSDLGLCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHSQWNKPSKPKTNMKHMAGAAAAGAVVGGLGGYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNITIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRGSSMVFFSSPPVILLISFLIFLIVG"]
    print(protein_bert_scores(ref, seqs))
        

if __name__ == "__main__":
    main()