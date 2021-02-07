#include "utils/unittest/Test.h"
#include "utils/multithread/MultiThread.h"
#include "utils/log/log.h"

TEST_MODULE_INIT(multi_thread_test);
using namespace BGIQD::LOG;
using namespace BGIQD::MultiThread;
TEST( multi_thread_run)
{
    timepoint start = timepoint::now();
    auto full_run = [] ()
    {
        int res = 0;
        for( int i = 0 ; i < 100000 ; i ++ )
            for( int j = 0 ; j < 10000 ; j++ )
                res +=j;
        return res;
    };

    BGIQD::MultiThread::MultiThread t_jobs;
    t_jobs.Start(3);
    t_jobs.AddJob( full_run);
    t_jobs.AddJob( full_run);
    t_jobs.AddJob( full_run);
    t_jobs.End();
    t_jobs.WaitingStop();
    timepoint end = timepoint::now();
    (test_logger)<<lstart()<<start.to_string()<<lend();
    (test_logger)<<lstart()<<end.to_string()<<lend();
    (test_logger)<<lstart()<<(end-start).to_string()<<lend();
}
