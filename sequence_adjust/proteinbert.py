

from protein_bert.proteinbert import load_pretrained_model
from protein_bert.proteinbert.conv_and_global_attention_model import get_model_with_hidden_layers_as_outputs

pretrained_model_generator, input_encoder = load_pretrained_model()

seqs = ["MQAQA"]
seq_len = len(seqs[0]) + 2
model = get_model_with_hidden_layers_as_outputs(pretrained_model_generator.create_model(seq_len))
encoded_x = input_encoder.encode_X(seqs, seq_len)
local_representations, global_representations = model.predict(encoded_x, batch_size=1)
# ... use these as features for other tasks, based on local_representations, global_representationss

print(global_representations.shape)
print(local_representations.shape)
