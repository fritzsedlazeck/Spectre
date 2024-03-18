FROM python:3.11.4-slim-buster

# Metadata
LABEL VERSION="0.2.0.alpha"
LABEL NAME="Spectre:${VERSION}"
LABEL AUTHOR="Philippe Sanio"



# Path: /app
WORKDIR /app/Spectre

# pip install requirements.txt
COPY requirements.txt .

# pip install
RUN pip install --no-cache-dir --upgrade pip \
  && pip install --no-cache-dir -r requirements.txt

RUN pip install -U setuptools wheel build

# copy Spectre directory
COPY . .

RUN python -m build

# Install what ever the dist folder contains and ends with .tar.gz
RUN pip install --no-cache-dir dist/*.tar.gz
RUN pip install --no-cache-dir .

# Start the application
ENTRYPOINT ["spectre"]
CMD ["Spectre"]