FROM sagemath/sagemath:9.7
USER root
RUN apt-get update && \
    apt-get install -y git && \
    apt-get clean
WORKDIR /app