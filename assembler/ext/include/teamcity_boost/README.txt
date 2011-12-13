Boost listener for TeamCity
---------------------------

To report your tests result to TeamCity server
just include teamcity_messages.* teamcity_boost.cpp
to your project.

That code will register global fixture
( http://www.boost.org/doc/libs/1_38_0/libs/test/doc/html/utf/user-guide/fixture/global.html )
to replace output formatter if run under TeamCity.

If you have tests with test parameters, see PARAM_TEST_CASES.txt for quick solution.

Technical details
-----------------

Reporting implemented as writing TeamCity service messages to stdout.

See
http://www.jetbrains.net/confluence/display/TCD3/Build+Script+Interaction+with+TeamCity
for more details.

Contact information
-------------------

See http://www.jetbrains.com/support/teamcity

License
-------

Apache, version 2.0
http://www.apache.org/licenses/LICENSE-2.0
