from flask import Flask
from polycip import Extraction

app = Flask(__name__)

@app.route('/')
def hello_world():
    example = 'FAMKOV.search1.mol2'
    response = Extraction.refine(example)
    return response

