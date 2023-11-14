ARG REG=100225593120.dkr.ecr.us-east-1.amazonaws.com
ARG DOCKER_PULL_TAG=stage

FROM ${REG}/agr_base_linux_env:${DOCKER_PULL_TAG}

WORKDIR /usr/src/app

ADD requirements.txt .

RUN . /root/venv/bin/activate

RUN pip install -r requirements.txt

RUN mkdir tmp

ADD . .

CMD ["python", "-u", "src/app.py", "--all-filetypes", "--upload"]
