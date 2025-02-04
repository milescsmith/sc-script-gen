FROM index.docker.io/python:3.12-slim

ENV DEBIAN_FRONTEND=noninteractive

COPY ./dist/sc_script_gen-2.6.0-py3-none-any.whl /opt/
RUN pip install /opt/sc_script_gen-2.6.0-py3-none-any.whl

ENTRYPOINT [ "script_gen" ]