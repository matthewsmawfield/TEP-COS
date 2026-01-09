#!/bin/bash
# =============================================================================
# IndexNow Submission Script for TEP-COS
# =============================================================================
# This script notifies search engines (Bing, Yandex, etc.) when the paper is updated.
# Run this after pushing changes to trigger immediate re-crawling.
#
# Usage: ./submit_indexnow.sh
# =============================================================================

# Configuration
HOST="matthewsmawfield.github.io"
KEY="29c6507763d2303d801cc8ed89d39f88"
KEY_LOCATION="https://mlsmawfield.com/tep/cos/${KEY}.txt"

# URLs to index
URL_TO_INDEX="https://mlsmawfield.com/tep/cos/"
# Note: Update PDF filename when available
PDF_TO_INDEX="https://mlsmawfield.com/tep/cos/public/docs/Smawfield_2026_TEP-COS_v0.1_Caracas.pdf"

# JSON Payload
JSON_PAYLOAD=$(cat <<EOF
{
  "host": "$HOST",
  "key": "$KEY",
  "keyLocation": "$KEY_LOCATION",
  "urlList": [
    "$URL_TO_INDEX",
    "$PDF_TO_INDEX"
  ]
}
EOF
)

echo "=============================================="
echo "IndexNow Submission for TEP-COS"
echo "=============================================="
echo ""
echo "Host: $HOST"
echo "Key Location: $KEY_LOCATION"
echo "URLs to index:"
echo "  - $URL_TO_INDEX"
echo "  - $PDF_TO_INDEX"
echo ""
echo "Submitting to IndexNow API..."
echo ""

# Submit to IndexNow (shared across Bing, Yandex, and other participating engines)
curl -s -X POST "https://api.indexnow.org/indexnow" \
     -H "Content-Type: application/json; charset=utf-8" \
     -d "$JSON_PAYLOAD" \
     -w "\nHTTP Status: %{http_code}\n"

echo ""
echo "=============================================="
echo "Submission complete!"
echo ""
echo "HTTP 200 = Success (URLs received by search engines)"
echo "HTTP 202 = Accepted (URLs queued for processing)"
echo "HTTP 400 = Bad Request (check JSON payload)"
echo "HTTP 403 = Forbidden (key validation failed)"
echo "HTTP 422 = Unprocessable (URLs don't match key location)"
echo "=============================================="
