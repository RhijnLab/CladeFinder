FROM mambaorg/micromamba:1.5.10
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes
WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml /app/
COPY --chown=$MAMBA_USER:$MAMBA_USER src /app/src
RUN micromamba run -n base python -m pip install --no-deps --no-cache-dir .
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "cladefinder"]
