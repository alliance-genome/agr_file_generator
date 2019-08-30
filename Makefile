build: pull
	docker build -t agrdocker/agr_file_generator:latest .

buildenv: build

pull:
	docker pull agrdocker/agr_base_linux_env:latest

run: build
	docker-compose up agr_file_generator

test: run
	docker-compose up agr_file_generator_test
