process SETUP {
    publishDir params.outdir, mode: 'symlink'
    conda "$projectDir/environment.yml"

    output:
    path "versions.json", emit: versions
    path "sawfish", emit: sawfish
    path "delly", emit: delly

    script:
    """
    wget https://github.com/PacificBiosciences/sawfish/releases/download/v0.12.9/sawfish-v0.12.9-x86_64-unknown-linux-gnu.tar.gz
    tar -xvf sawfish-v0.12.9-x86_64-unknown-linux-gnu.tar.gz
    mv sawfish-v0.12.9-x86_64-unknown-linux-gnu/bin/sawfish .

    wget https://github.com/dellytools/delly/releases/download/v1.3.3/delly_v1.3.3_linux_x86_64bit
    mv delly_v1.3.3_linux_x86_64bit delly && chmod +x delly

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
