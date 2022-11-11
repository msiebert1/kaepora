"""test Flask with this"""

from flask import Flask, render_template, request
from flask_sqlalchemy import SQLAlchemy
from flask_bootstrap import Bootstrap
from sqlalchemy.sql import text
import sys,os

sys.path.insert(1, '../src')
import composite
import kaepora as kpora

app = Flask(__name__)
Bootstrap(app)

db_name = '../data/kaepora_v1.2.db'


print(f'os.path.abspath(db_name): {str(os.path.abspath(db_name))}')
app.config['SQLALCHEMY_DATABASE_URI'] = "sqlite:///" + os.path.abspath(db_name)
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = True
db = SQLAlchemy(app)


# each table in the database needs a class to be created for it
# db.Model is required - don't change it
# identify all columns by name and data type
class Spectrum(db.Model):
    __tablename__ = 'spectra'
    filename = db.Column(db.String, primary_key=True)
    sn = db.Column(db.String)

class Event(db.Model):
    __tablename__ = 'events'
    sn = db.Column(db.String, primary_key=True)
    redshift = db.Column(db.String)


# @app.route('/')
# def index():
#     try:
#         # spectra = Spectrum.query.all()
#         # spec_text = '<ul>'
#         # for i, spec in enumerate(spectra):
#         #     # spec_text += '<li>' + spec.SN + '</li>'
#         #     spec_text += '<li>' + str(i) + '  ' + spec.SN + ', ' + spec.filename + '</li>'
#         # spec_text += '</ul>'
#         # return spec_text

#         sne = Event.query.all()
#         sn_text = '<ul>'
#         for i, supernova in enumerate(sne):
#             # spec_text += '<li>' + spec.SN + '</li>'
#             sn_text += '<li>' + str(i) + '  ' + supernova.sn + ', ' + str(supernova.redshift) + '</li>'
#         sn_text += '</ul>'
#         return sn_text

#     except Exception as e:
#         # e holds description of the error
#         error_text = "<p>The error:<br>" + str(e) + "</p>"
#         hed = '<h1>Something is broken.</h1>'
#         return hed + error_text


@app.route('/foobar')
def foobar():
    return '<h1>Hi there, foobar!</h1>'


@app.route('/')
def inventory():
    spectra = Spectrum.query.all()
    return render_template('list.html', spectra=spectra)

@app.route('/form')
def form():
    return render_template('form.html')
 
@app.route('/data/', methods = ['POST', 'GET'])
def data():
    if request.method == 'GET':
        return f"The URL /data is accessed directly. Try going to '/form' to submit form"
    if request.method == 'POST':
        form_data = request.form
        phase_low, phase_up = form_data['phase'].split()
        dm15_low, dm15_up = form_data['dm15'].split()

        query = "SELECT * from Spectra inner join Events ON Spectra.SN = Events.SN where" \
                 " phase >= {phase_low} and phase <= {phase_up} and ((dm15_source between {dm15_low}" \
                 " and {dm15_up}) or (dm15_from_fits between {dm15_low} and {dm15_up}))".format(phase_low=phase_low, phase_up=phase_up, dm15_low=dm15_low, dm15_up=dm15_up)
        print (query)
        spec_array = kpora.grab(query, verbose=False, db_file = db_name)
        spec_text = '<ul>'
        for spec in spec_array:
            # spec_text += '<li>' + spec.SN + '</li>'
            spec_text += '<li>' + spec.name + ', ' + spec.filename + '</li>'
        spec_text += '</ul>'
        return spec_text
        # return render_template('data.html',form_data = form_data)

# @app.route('/')
# def index():
#     try:
#         example_query = ["SELECT * from Spectra inner join Events ON Spectra.SN = Events.SN where" \
#                          " phase >= -6 and phase <= -4 and ((dm15_source between .5 and 1.8) or (dm15_from_fits between .5 and 1.8))"]

#         # spectra = Spectrum.query.filter_by(SN='2011fe').order_by(Spectrum.SN).all()
#         spec_array = kpora.grab(example_query[0], verbose=False, db_file = db_name)
#         spec_text = '<ul>'
#         for spec in spec_array:
#             # spec_text += '<li>' + spec.SN + '</li>'
#             spec_text += '<li>' + spec.name + ', ' + spec.filename + '</li>'
#         spec_text += '</ul>'
#         return spec_text
#     except Exception as e:
#         # e holds description of the error
#         error_text = "<p>The error:<br>" + str(e) + "</p>"
#         hed = '<h1>Something is broken.</h1>'
#         return hed + error_text


if __name__ == '__main__':
    app.run(debug=True)





