if ~exist('../../build/bin', 'dir')
	error('Test binaries are missing. Use CMake to build the tests.')
end

test__rk4_solver__sine
test__rk4_solver__motor