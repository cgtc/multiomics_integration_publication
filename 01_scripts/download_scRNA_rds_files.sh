#!/bin/bash

# Script to download raw fastq.gz files to the local computer for testing the Nextflow pipeline

echo "Script started..."

# SSH key for connecting to the AWS EC2 instance
EC2_SSH_KEY="C:\Users\joddy\OneDrive - Cell and Gene Therapy Catapult\aws\awsKeyMFastenrath.pem"

# AWS EC2 instance details
EC2_USER="ubuntu"
EC2_HOST="ec2-3-11-245-196.eu-west-2.compute.amazonaws.com"

# Base directory on the EC2 instance
BASE_DIR="/efs/plasticell/integration/plasticell_multiomics_integration/00_data/raw/unprocessed/scRNA_processing_files"

# Output directory on your local machine
OUTPUT_DIR="/c/Users/joddy/repos/plasticell_multiomics_integration/00_data/raw/unprocessed/scRNA_processing_files"

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Get the list of rds files from EC2 instance
echo "Getting list of rds files from the EC2 instance..."
FILES=$(ssh -i "$EC2_SSH_KEY" "$EC2_USER@$EC2_HOST" "find $BASE_DIR -type f -name '*.rds'")

# Download each file using scp with progress
for FILE in $FILES
do
    # Extract the file name
    FILE_NAME=$(basename "$FILE")
    
    # Display progress for each file
    echo "Downloading $FILE_NAME..."

    # Use SCP to download the file and show progress
    scp -i "$EC2_SSH_KEY" "$EC2_USER@$EC2_HOST:$FILE" "$OUTPUT_DIR/$FILE_NAME" &
    
    # Wait for the background SCP process to complete
    wait $!
    
    # Inform the user that the file has been downloaded
    echo "$FILE_NAME downloaded successfully."
done

echo "All files downloaded. Files are saved in $OUTPUT_DIR."