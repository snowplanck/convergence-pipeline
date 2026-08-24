# ML Structure Annotation Stack (B3 — GPU)
# ProstT5, ESMFold, DeepFRI, CLEAN each have distinct CUDA/PyTorch requirements.
# Build separate images per tool to isolate dependency stacks.

FROM pytorch/pytorch:2.1.0-cuda12.1-cudnn8-runtime

RUN pip install transformers accelerate sentencepiece biopython

# ProstT5
RUN pip install torch-scatter -f https://data.pyg.org/whl/torch-2.1.0+cu121.html

# DeepFRI (graph convolutional network)
RUN pip install dgl -f https://data.dgl.ai/wheels/cu121/repo.html
RUN pip install deepfri

# ESMFold
RUN pip install fair-esm

# CLEAN
RUN pip install clean-nn

WORKDIR /app
ENTRYPOINT ["python3"]
