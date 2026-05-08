# Use a lightweight Python 3.10 image
FROM python:3.10-slim

# Set the working directory inside the container
WORKDIR /app

# Install system-level dependencies for Matplotlib and Biopython
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy the requirements file first to leverage Docker cache
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of your application code
COPY . .

# Render uses port 10000 by default for Web Services
EXPOSE 10000

# Run the application
# We use uvicorn logic through app.py's demo.launch()
CMD ["python", "app.py"]