# MitoM mitogenome annotation workflow (one conda environment, env/mitotools.yml, + MITOS2's refseq89m) for any
# system with Docker (Windows, macOS incl. Apple Silicon via emulation, Linux). Built by install.sh in a clean Ubuntu.
#   docker build -t mitom .
#   docker run --rm -v "$PWD":/data mitom mitom run -i /data/genomes.gb -o /data/out
FROM --platform=linux/amd64 ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive LANG=C.UTF-8
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates curl bzip2 procps \
        && rm -rf /var/lib/apt/lists/*
COPY . /opt/mito
RUN bash /opt/mito/install.sh --refdir /opt/mitos_refs && /root/miniforge3/bin/mamba clean -a -y
# the mitom command (install.sh) runs on the environment's Python
ENV PATH=/root/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
# build-time check: the command starts and every required program and the reference data are found
RUN mitom --version && mitom run -h >/dev/null && mitom doctor
WORKDIR /data
CMD ["mitom", "doctor"]
