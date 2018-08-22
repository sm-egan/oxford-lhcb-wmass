#include <iostream>
#include <string>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int main( int argc, const char** argv )
{
  std::string a_string_option;
  int an_int_option;
  
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("int_option", po::value(&an_int_option)->default_value(1))
    ("string_option", po::value(&a_string_option)->default_value("hello"))
    ;
  
  po::variables_map vm;
  po::store( po::parse_command_line( argc, argv, desc ), vm );
  if ( vm.count( "help" ) ) {
    std::cout << desc << std::endl;
    return 1;
  }
  po::notify( vm );
  
  return 0;
}
