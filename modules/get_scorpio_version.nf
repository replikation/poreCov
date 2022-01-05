process get_scorpio_version {
    label 'pangolin'
    cpus = 1
    container = params.pangolindocker
  output:
    tuple env(SCORPIO_VER), env(SCORPIO_CONSTELLATIONS_VER)
  shell:
    '''
    SCORPIO_VER=$(scorpio --version)
    SCORPIO_CONSTELLATIONS_VER=$(scorpio -cv)
    '''
    stub:
    """
    SCORPIO_VER='scorpio 0.3.14'
    SCORPIO_CONSTELLATIONS_VER='constellations v0.0.24'
    """
  }
