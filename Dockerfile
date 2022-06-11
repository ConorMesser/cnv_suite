FROM gcr.io/broad-getzlab-workflows/base_image:v0.0.5

WORKDIR /build

# install cnv_suite
COPY setup.py .
COPY cnv_suite ./cnv_suite
RUN pip install .

WORKDIR /app
ENV PATH=$PATH:/app
