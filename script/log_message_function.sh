#!/bin/bash

# --- Color Definitions ---
# \033 or \e starts the escape sequence
# [0m resets all formatting
RED='\033[0;31m'
YELLOW='\033[0;33m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color - Reset all attributes

# Bold text for the log type
BOLD_RED='\033[1;31m'
BOLD_YELLOW='\033[1;33m'
BOLD_GREEN='\033[1;32m'
BOLD_BLUE='\033[1;34m'

# --- Logging Functions ---

log_info() {
  echo "#############"
  echo -e "${BOLD_BLUE}[INFO]${NC} ${GREEN}$1${NC}"
  echo "#############"

}

log_warning() {
  echo "#############"
  echo -e "${BOLD_YELLOW}[WARNING]${NC} ${YELLOW}$1${NC}"
  echo "#############"

}

log_error() {
  # Directs to stderr (standard error)
  echo "#############"
  echo -e "${BOLD_RED}[ERROR]${NC} ${RED}$1${NC}" >&2
  echo "#############"

}
