

FROM pmpnn_env:latest


COPY ./ProteinMPNN/ /workspace/pmpnn/

RUN chmod -R 777 /workspace/pmpnn/




