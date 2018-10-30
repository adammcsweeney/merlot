# module for configuring settings for WORHP
# for more details refer to the users guide: https://worhp.de/latest/download/user_manual.pdf

# TODO ==== perform study to find ideal operating values, i.e. comprimising run time and likelihood of optimal trajectory being found

def worhp_setup(uda):

    # WORHP setup

    # TolFeas
    # Controls the tolerance to be imposed on each constraint to be regarded as feasible (default: 1e-6)
    uda.set_numeric_option('TolFeas', 1e-6)

    # TolOpti
    # Controls the tolerance to be imposed on the (scaled) KKT conditions for a point to be
    # regarded as optimal (default 1e-6)
    uda.set_numeric_option('TolOpti', 1e-6)

    # TolComp
    # Controls the complementarity tolerance to be imposed on the multipliers (default 1e-6)
    uda.set_numeric_option('TolComp', 1e-6)

    # # AcceptTolFeas
    # # Controls the acceptable feasibility tolerance.  This value is computed during initialisation depending on TolFeas,
    # # but may be overridden by the user.  Setting it to the same value as its strict counterpart disables the check for
    # # acceptable feasibility i.e. must be >= TolFeas (default 1e-3)
    # uda.set_numeric_option('AcceptTolFeas', 1e-3)
    #
    # # AcceptTolOpti
    # # Controls the acceptable optimality tolerance.  This value is computed during initialisation depending on TolOpti,
    # # but may be overridden by the user.  Setting it to the same value as its strict counterpart disables the check for
    # # acceptable optimality >= TolOpti (default 1e-3)
    # uda.set_numeric_option('AcceptTolOpti', 1e-3)

    # MaxIter
    # Imposes an upper limit on the number of major solver iterations.  If the limit is reached, Worhp will check for an
    # acceptable point and terminate (default 10000)
    uda.set_integer_option('MaxIter', 500)

    return uda