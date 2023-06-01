#!/bin/bash

set -euo pipefail
# set -euxo pipefail

FAEST_128S="faest_128s"
FAEST_128F="faest_128f"
FAEST_192S="faest_192s"
FAEST_192F="faest_192f"
FAEST_256S="faest_256s"
FAEST_256F="faest_256f"
FAEST_EM_128S="faest_em_128s"
FAEST_EM_128F="faest_em_128f"
FAEST_EM_192S="faest_em_192s"
FAEST_EM_192F="faest_em_192f"
FAEST_EM_256S="faest_em_256s"
FAEST_EM_256F="faest_em_256f"

ALL="FAEST_128S FAEST_128F FAEST_192S FAEST_192F FAEST_256S FAEST_256F FAEST_EM_128S FAEST_EM_128F FAEST_EM_192S FAEST_EM_192F FAEST_EM_256S FAEST_EM_256F"

BUILD_DIR="build_release"

clean () {
    make clean
    (cd ${BUILD_DIR} && ninja clean)
}

run_make () {
    name="$1"
    setting_id="$2"
    echo "# Building: ${name} (${setting_id})" >&2
    (cd ${BUILD_DIR} && ninja "${setting_id}/${setting_id}_api_test") >&2
    (cd ${BUILD_DIR} && ninja "${setting_id}/${setting_id}_bench_c2") >&2
}

run_test () {
    name="$1"
    setting_id="$2"
    echo "# Testing: ${name} (${setting_id})" >&2
    "./${BUILD_DIR}/${setting_id}/${setting_id}_api_test" >&2
}

convert_to_us () {
    input="$1"
    if ! echo "${input}" | grep -Pq '^\d+(\.\d+)? (u|m|n)?s$'; then
        echo "unexpected format: '${input}'" >&2
        exit 1
    fi
    num="`echo "${input}" | cut -d ' ' -f 1`"
    unit="`echo "${input}" | cut -d ' ' -f 2`"
    if [ "$unit" = "ns" ]; then
        echo "scale=10; ${num} / 1000"  | bc
    elif [ "$unit" = "us" ]; then
        echo "${num}"
    elif [ "$unit" = "ms" ]; then
        echo "scale=10; ${num} * 1000"  | bc
    elif [ "$unit" = "s" ]; then
        echo "scale=10; ${num} * 1000 * 1000"  | bc
    fi
}

parse_output () {
    output="$1"
    line1="`echo "${output}" | sed -e '1q;d' | tr -s ' '`"
    line2="`echo "${output}" | sed -e '2q;d' | tr -s ' '`"
    line3="`echo "${output}" | sed -e '3q;d' | tr -s ' '`"
    samples="`echo ${line1} | cut -d ' ' -f 2`"
    iterations="`echo ${line1} | cut -d ' ' -f 3`"
    estimated="`echo ${line1} | cut -d ' ' -f 4-5`"
    mean="`echo ${line2} | cut -d ' ' -f 1-2`"
    low_mean="`echo ${line2} | cut -d ' ' -f 3-4`"
    high_mean="`echo ${line2} | cut -d ' ' -f 5-6`"
    std_dev="`echo ${line3} | cut -d ' ' -f 1-2`"
    low_std_dev="`echo ${line3} | cut -d ' ' -f 3-4`"
    high_std_dev="`echo ${line3} | cut -d ' ' -f 5-6`"
    json="`jq -n \
        --argjson "samples" "${samples}" \
        --argjson "iterations" "${iterations}" \
        --argjson "estimated" "$(convert_to_us "${estimated}")" \
        --argjson "mean_us" "$(convert_to_us "${mean}")" \
        --argjson "low_mean_us" "$(convert_to_us "${low_mean}")" \
        --argjson "high_mean_us" "$(convert_to_us "${high_mean}")" \
        --argjson "std_dev_us" "$(convert_to_us "${std_dev}")" \
        --argjson "low_std_dev_us" "$(convert_to_us "${low_std_dev}")" \
        --argjson "high_std_dev_us" "$(convert_to_us "${high_std_dev}")" \
        '$ARGS.named' \
        `"
    echo "${json}"
}

gather_metadata () {
    hostname="`uname -n`"
    user="`whoami`"
    timestamp="`date --iso-8601=ns`"
    git_commit="`git rev-parse --verify HEAD`"
    json="`jq -c -n \
        --arg "hostname" "${hostname}" \
        --arg "user" "${user}" \
        --arg "timestamp" "${timestamp}" \
        --arg "git_commit" "${timestamp}" \
        '$ARGS.named' \
        `"
    echo "${json}"
}

run_bench () {
    name="$1"
    setting_id="$2"
    echo "# Benching: ${name} (${setting_id})" >&2

    meta="`gather_metadata`"

    bench_out="`"./${BUILD_DIR}/${setting_id}/${setting_id}_bench_c2" '[bench]'`"

    bench_keygen_out="`echo "${bench_out}" | grep -P '^keygen' -A 2`"
    bench_sign_out="`echo "${bench_out}" | grep -P '^sign' -A 2`"
    bench_verify_out="`echo "${bench_out}" | grep -P '^verify' -A 2`"

    sk_size="`grep CRYPTO_SECRETKEYBYTES < "./${BUILD_DIR}/${setting_id}/api.h" | cut -d ' ' -f 3- | bc`"
    pk_size="`grep CRYPTO_PUBLICKEYBYTES < "./${BUILD_DIR}/${setting_id}/api.h" | cut -d ' ' -f 3`"
    sig_size="`grep CRYPTO_BYTES < "./${BUILD_DIR}/${setting_id}/api.h" | cut -d ' ' -f 3`"
    keygen_results="`parse_output "${bench_keygen_out}"`"
    sign_results="`parse_output "${bench_sign_out}"`"
    verify_results="`parse_output "${bench_verify_out}"`"

    json="`jq -c -n \
        --arg "implementation" "ref" \
        --arg "variant" "${name}" \
        --arg "setting_id" "${setting_id}" \
        --argjson "sig_size_bytes" "${sig_size}" \
        --argjson "sk_size_bytes" "${sk_size}" \
        --argjson "pk_size_bytes" "${pk_size}" \
        --argjson "keygen" "${keygen_results}" \
        --argjson "sign" "${sign_results}" \
        --argjson "verify" "${verify_results}" \
        --argjson "meta" "${meta}" \
        '$ARGS.named' \
        `"
    echo "${json}"
}

bench_spec_variants () {
    for faest_variant in $ALL; do
        declare -n setting_id="${faest_variant}"
        run_make "${faest_variant}" "${setting_id}"
        run_test "${faest_variant}" "${setting_id}"
        run_bench "${faest_variant}" "${setting_id}"
    done
}

bench_spec_variants
