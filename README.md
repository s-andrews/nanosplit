# NanoSplit

NanoSplit is a Flask web application for splitting mixed nanopore FASTQ reads by barcode. It guides a user through upload, read-length filtering, barcode discovery, final barcode or barcode-pair quantitation, and download of split FASTQ files.

## Functionality

NanoSplit provides a browser-based workflow:

1. Upload a `.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz` file.
2. Choose barcode length and whether reads should be matched at one end or both ends.
3. Review the read-length distribution and submit a length filter.
4. Review barcode abundance at the start of reads.
   - Barcode sequences are plotted as horizontal bars.
   - Selected barcodes can be controlled directly from the plot.
   - A selected barcode also represents its reverse complement.
5. Review final barcode or barcode-pair abundance.
   - Paired barcode results merge opposite orientations into stacked forward/reverse bars.
   - Barcode pairs can include the same barcode at both ends.
6. Download split FASTQ files individually or as a ZIP archive.

Job files are stored under the local `data/` directory. Output filenames are based on the uploaded filename, with barcode names and sequences appended where available.

## Local Development

Install Python 3.10 or newer, then create a virtual environment:

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Run the development server:

```bash
python app.py
```

Open `http://127.0.0.1:5000/`.

## Linux Server Installation

These commands assume the application will live at `/opt/nanosplit`, run as user `www-data`, and listen locally on port `8050`.

```bash
sudo apt update
sudo apt install python3 python3-venv apache2
sudo mkdir -p /opt/nanosplit
sudo chown www-data:www-data /opt/nanosplit
```

Copy the NanoSplit repository into `/opt/nanosplit`, then install dependencies:

```bash
cd /opt/nanosplit
sudo -u www-data python3 -m venv venv
sudo -u www-data ./venv/bin/pip install -r requirements.txt
sudo -u www-data ./venv/bin/pip install gunicorn
sudo -u www-data mkdir -p data
```

## Run With systemd

Create `/etc/systemd/system/nanosplit.service`:

```ini
[Unit]
Description=NanoSplit web application
After=network.target

[Service]
User=www-data
Group=www-data
WorkingDirectory=/opt/nanosplit
Environment=PYTHONUNBUFFERED=1
ExecStart=/opt/nanosplit/venv/bin/gunicorn --workers 2 --bind 127.0.0.1:8050 app:app
Restart=always

[Install]
WantedBy=multi-user.target
```

Start the service:

```bash
sudo systemctl daemon-reload
sudo systemctl enable --now nanosplit
sudo systemctl status nanosplit
```

The site should now be listening only on the server itself at `http://127.0.0.1:8050/`.

## Apache Reverse Proxy at `/nanosplit/`

Enable Apache proxy modules:

```bash
sudo a2enmod proxy proxy_http headers
sudo systemctl reload apache2
```

Add this to the relevant Apache virtual host, for example `/etc/apache2/sites-available/example.com.conf`:

```apache
<VirtualHost *:443>
    ServerName example.com

    # Existing TLS configuration goes here.

    ProxyPreserveHost On
    RequestHeader set X-Forwarded-Prefix "/nanosplit"
    RequestHeader set X-Forwarded-Proto "https"

    ProxyPass        /nanosplit/ http://127.0.0.1:8050/
    ProxyPassReverse /nanosplit/ http://127.0.0.1:8050/
</VirtualHost>
```

Reload Apache:

```bash
sudo apachectl configtest
sudo systemctl reload apache2
```

NanoSplit should then be available at:

```text
https://example.com/nanosplit/
```

To use a different local port, change both the `--bind 127.0.0.1:8050` value in the systemd service and the Apache `ProxyPass` / `ProxyPassReverse` target.

## Operational Notes

- Uploaded and generated files are kept in `/opt/nanosplit/data`.
- Ensure that the service user can write to `data/`.
- Apache handles public HTTPS traffic; Gunicorn should remain bound to `127.0.0.1`.
- For large FASTQ files, configure Apache upload limits and server disk capacity appropriately.
- The development server in `python app.py` is not intended for production hosting.
