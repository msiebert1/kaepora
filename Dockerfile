# syntax=docker/dockerfile:1

FROM python:3.10.6

ENV FLASK_APP=./flask_app/kaepora_app.py

COPY ./requirements.txt /requirements.txt
RUN pip install -r requirements.txt

EXPOSE 5000

COPY . .

CMD [ "python3", "-m" , "flask", "run", "--host=0.0.0.0"]