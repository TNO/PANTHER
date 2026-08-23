classdef TestMultiFaultAnalyzer < matlab.unittest.TestCase
    % TestMultiFaultAnalyzer Functional tests for MultiFaultAnalyzer controller.
    % Tests verify that running multiple faults via MultiFaultAnalyzer
    % produces the same results as running each FaultAnalyzer individually.

    methods (Test)

        function test_sequential_run_matches_single(testCase)
            % A MultiFaultAnalyzer with 2 identical default faults run
            % sequentially should produce the same sne/tau as a single
            % FaultAnalyzer run.
            ref = FaultAnalyzer();
            ref = ref.run();

            mfa = MultiFaultAnalyzer();
            mfa = mfa.initialize(2);
            mfa.parallel = 0;
            mfa.printStatusOutput = false;
            mfa = mfa.run();

            i_mid = ceil(length(ref.y) / 2);
            for k = 1 : mfa.nFaults
                actual_sne = mfa.faults(k).faultResults.sne(i_mid, end);
                testCase.verifyEqual(actual_sne, ref.faultResults.sne(i_mid, end), ...
                    'RelTol', 1e-6, ...
                    sprintf('sne mismatch for fault %d', k));
                actual_tau = mfa.faults(k).faultResults.tau(i_mid, end);
                testCase.verifyEqual(actual_tau, ref.faultResults.tau(i_mid, end), ...
                    'RelTol', 1e-6, ...
                    sprintf('tau mismatch for fault %d', k));
            end
        end

        function test_parallel_matches_sequential(testCase)
            % Parallel run should give identical results to sequential run.
            mfa_seq = MultiFaultAnalyzer();
            mfa_seq = mfa_seq.initialize(3);
            mfa_seq.parallel = 0;
            mfa_seq.printStatusOutput = false;
            mfa_seq = mfa_seq.run();

            mfa_par = MultiFaultAnalyzer();
            mfa_par = mfa_par.initialize(3);
            mfa_par.parallel = 1;
            mfa_par.printStatusOutput = false;
            mfa_par = mfa_par.run();

            i_mid = ceil(length(mfa_seq.faults(1).y) / 2);
            for k = 1 : mfa_seq.nFaults
                testCase.verifyEqual( ...
                    mfa_par.faults(k).faultResults.sne(i_mid, end), ...
                    mfa_seq.faults(k).faultResults.sne(i_mid, end), ...
                    'RelTol', 1e-6, ...
                    sprintf('sne parallel/sequential mismatch for fault %d', k));
            end
        end

        function test_run_done_flag(testCase)
            % runDone should be true after run().
            mfa = MultiFaultAnalyzer();
            mfa = mfa.initialize(2);
            mfa.parallel = 0;
            mfa.printStatusOutput = false;
            mfa = mfa.run();
            testCase.verifyTrue(logical(mfa.runDone));
        end

        function test_fault_summary_populated(testCase)
            % faultSummary should be a non-empty table after run().
            mfa = MultiFaultAnalyzer();
            mfa = mfa.initialize(2);
            mfa.setInputParameter('dip',[60,70]);
            mfa.parallel = 0;
            mfa.printStatusOutput = false;
            mfa = mfa.run();
            testCase.verifyNotEmpty(mfa.faultSummary);
            testCase.verifyTrue(istable(mfa.faultSummary));
            testCase.verifyEqual(height(mfa.faultSummary), mfa.nFaults);
        end

        function test_get_result_for_multiple_faults(testCase)
            mfa = MultiFaultAnalyzer();
            mfa = mfa.initialize(2);
            mfa.parallel = 0;
            mfa.printStatusOutput = false;
            mfa = mfa.run();

            results = mfa.getResult('sne');
            testCase.verifyClass(results, 'cell');
            testCase.verifySize(results, [2, 1]);
            testCase.verifyEqual(results{1}, mfa.faults(1).faultResults.sne);
            testCase.verifyEqual(mfa.getResult('tau', 2), mfa.faults(2).faultResults.tau);
            testCase.verifyEqual(mfa.getResultAtLoadStep('sne', 1, 2), mfa.faults(2).faultResults.sne(:, 1));
            testCase.verifyEqual(mfa.getResultAtY('sne', mfa.faults(1).y(1), 1), mfa.faults(1).faultResults.sne(1, :));
        end

        function test_print_status_every_n(testCase)
            % printStatusEveryNFaults should not cause errors.
            mfa = MultiFaultAnalyzer();
            mfa = mfa.initialize(5);
            mfa.parallel = 0;
            mfa.printStatusOutput = true;
            mfa.printStatusEveryNFaults = 2;
            mfa = mfa.run();
            testCase.verifyEqual(mfa.nFaults, 5);
        end

    end
end


