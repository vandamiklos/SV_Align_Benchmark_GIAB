/*process SETUP {
    publishDir params.outdir, mode: 'symlink'
    conda "$projectDir/environment.yml"

    output:
    path "versions.json", emit: versions
    path "sawfish", emit: sawfish
    path "delly", emit: delly

    script:
    """
#    wget https://github.com/PacificBiosciences/sawfish/releases/download/v0.12.9/sawfish-v0.12.9-x86_64-unknown-linux-gnu.tar.gz
#    tar -xvf sawfish-v0.12.9-x86_64-unknown-linux-gnu.tar.gz
#    mv sawfish-v0.12.9-x86_64-unknown-linux-gnu/bin/sawfish .

#    wget https://github.com/dellytools/delly/releases/download/v1.3.3/delly_v1.3.3_linux_x86_64bit
#    mv delly_v1.3.3_linux_x86_64bit delly && chmod +x delly

    dysgu_ver=\$(dysgu --version | cut -f3 -d ' ')
    sniffles_ver=\$(sniffles --version | cut -f3 -d ' ')
    cutesv_ver=\$(cuteSV --version | cut -f2 -d ' ')
    sawfish_ver=\$(./sawfish --version | cut -f2 -d ' ')

    cat <<- EOF > versions.json
    {
        "dysgu": "\$dysgu_ver",
        "sniffles": "\$sniffles_ver",
        "cutesv": "\$cutesv_ver",
        "sawfish": "\$sawfish_ver",
        "delly": "1.3.3"
    }
EOF
    """
}



// chatgpt
process SETUP {

    tag "setup"

    input:
    path ref

    output:
    path "versions.json", emit: versions

    script:
    """
    # Collect versions (robust)
    cat <<- EOF > versions.json
    {
        "samtools": "$(samtools --version | head -n1 | awk '{print $2}')",
        "sniffles": "$(sniffles --version 2>/dev/null | awk '{print $NF}')",
        "cutesv": "$(cuteSV --version 2>/dev/null | awk '{print $NF}')",
        "dysgu": "$(dysgu --version 2>/dev/null | awk '{print $NF}')",
        "delly": "$(delly 2>&1 | grep Version | awk '{print $2}')",
        "sawfish": "$(sawfish --version 2>/dev/null | awk '{print $NF}')"
    }
EOF
    """
} */