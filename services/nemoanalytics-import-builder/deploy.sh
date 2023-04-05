#!/usr/bin/env bash
set -euo pipefail

declare -r SCRIPT_DIR=$(cd $(dirname $0) && pwd)
declare -r CORS_CONFIG=${SCRIPT_DIR}/cors-config.json

# What we really want is to have a set of valid environment tokens
# with a "contains" check.
# As far as I can tell that isn't supported in bash.
# Instead we use the valid env strings as map keys with non-empty values.
# To check if a token is a valid env, we index into the map and assert
# that the result is non-empty.
declare -rA VALID_ENVS=([dev]=valid [prod]=valid)

function check_usage () {
  if [[ $# -ne 1 ]]; then
    2>&1 echo Error: Incorrect number of arguments given, expected 1 '(environment)' but got $#
    exit 1
  elif [[ -z "${VALID_ENVS[$1]-}" ]]; then
    2>&1 echo Error: Invalid environment "'$1'", valid values are: ${!VALID_ENVS[@]}
    exit 1
  fi
}

function secret_path () {
    echo "secret/dsde/monster/$1/biccn/nemoanalytics-import-builder/$2"
}

function set_cors () {
  gsutil cors set ${CORS_CONFIG} "gs://$1"
}

function deploy () {
    gcloud --project=$1 functions deploy \
        build-nemoanalytics-import \
        --runtime=nodejs6 \
        --trigger-http \
        --entry-point=buildNemoanalyticsImport \
        --set-env-vars TMP_BUCKET=$2 \
        --source=${SCRIPT_DIR}/function \
        --service-account=$3
}

main () {
    check_usage ${@}
    local -r env=$1

    local -r env_secrets=$(secret_path ${env} env)
    local -r account_secrets=$(secret_path ${env} service-account.json)

    local -r deploy_project=$(vault read -field=project ${env_secrets})
    local -r tmp_bucket=$(vault read -field=bucket ${env_secrets})
    local -r service_account_email=$(vault read -field=client_email ${account_secrets})

    set_cors ${tmp_bucket}
    deploy ${deploy_project} ${tmp_bucket} ${service_account_email}
}

main ${@}
