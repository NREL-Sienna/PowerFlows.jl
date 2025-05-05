@enum MUMPS_JOB begin
    INIT = -1
    CLEANUP = -2
    ANALYZE = 1
    FACTOR = 2
    SOLVE = 3
    SAVE = 7
    RESTORE = 8
    DELETE = -3
    FACTOR_CLEANUP = -4
end

function mumps_job!(mumps::Mumps, job::MUMPS_JOB)
    MUMPS.set_job!(mumps, Integer(job))
    MUMPS.invoke_mumps!(mumps)
    return
end
