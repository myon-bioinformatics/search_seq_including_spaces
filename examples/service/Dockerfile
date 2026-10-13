FROM python:3.13-slim
WORKDIR /app
COPY *.py ./
COPY protein ./protein
USER 65534:65534
EXPOSE 8765
CMD ["python", "-S", "sequence_service.py", "http", "--host", "0.0.0.0"]
