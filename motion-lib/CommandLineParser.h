// ====================================================================================================   
//                                                                                                        
//   SuPReMo: Surrogate Parameterised Respiratory Motion Model                                            
//            An implementation of the generalised motion modelling and image registration framework      
//                                                                                                        
//   Copyright (c) University College London (UCL). All rights reserved.                                  
//                                                                                                        
//   This software is distributed WITHOUT ANY WARRANTY; without even                                      
//   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR                                  
//   PURPOSE.                                                                                             
//                                                                                                        
//   See LICENSE.txt in the top level directory for details.                                              
//                                                                                                        
// ====================================================================================================   




#pragma once

#include <map>
#include <string>
#include <vector>
#include "Supremo.h"

struct CommandLineOption{
  unsigned int numberOfComponents; ///< Number of elements expected after this flag on the command line
  bool required;                   ///< Indicating if this parameters is essential
  std::string description;         ///< Description of the parameters
};


/** Function to split a string by a delimiter into a vector of strings. Items will only be added if they are not empty. 
 *  I.e. two commas will not add an empty element. 
 * 
 *  \param stringToSplitIn string which will be split
 *  \param delimiter Delimiter which will be used to split the string
 */
std::vector<std::string> splitStringbyDelimiter(const std::string & stringToSplitIn, const std::string & delimiter);


/** Class implementing a simple command line parser
 * Command line parser (heavily) adapted from 
 * http://stackoverflow.com/questions/865668/ddg#868894
 * credit to iain.
 */
class CommandLineParser{
public:

  /** Constructor.
   *  Saves each item of the command line in the tokens member variable and the executable name in executableName 
   *
   *  \param argc Number of input arguments.
   *  \param argv Pointer ot character arrays.
   */
  CommandLineParser( int &argc, char **argv, std::map<std::string, CommandLineOption> & allowedCommandLineOptionsIn );

  /** Get an option associated with a certain flag.
   *
   * Example
   * -------
   *  If the application was called like @verbatim application.exe -f fileName.nii.gz -dynamic 150 textFile.txt @endverbatim
   *  then @verbatim CommandLineParser.getCmdOption("-f")          @endverbatim returns "fileName.nii.gz",
   *  and  @verbatim CommandLineParser.getCmdOption("-dynamic", 0) @endverbatim returns "150",
   *  and  @verbatim CommandLineParser.getCmdOption("-dynamic", 1) @endverbatim returns "textFile.txt",
   *
   * \param  option String with the flag that is sought.
   * \param  skip   Defines how many tokens to skip after flag before selecting return value.
   * \return String with the specified option. Will be empty if option does not exist.
   */
  const std::string& getCmdOptionAsString( const std::string &option, int skip = 0 ) const;

  /** Get an option associated with a certain flag as an integer.
  *   For an example see getCmdOptionAsString.
  * \param  option String with the flag that is sought.
  * \param  skip   Defines how many tokens to skip after flag before selecting return value.
  * \return Integer with the specified option. Will exit if option does not exist as conversion is undefined in this case.
  */
  const int getCmdOptionAsInt( const std::string &option, int skip = 0 ) const;

  /** Get an option associated with a certain flag as a float.
  *   For an example see getCmdOptionAsString.
  * \param  option String with the flag that is sought.
  * \param  skip   Defines how many tokens to skip after flag before selecting return value.
  * \return Integer with the specified option. Will exit if option does not exist as conversion is undefined in this case.
  */
  const float getCmdOptionAsFloat( const std::string &option, int skip = 0 ) const;

  /**
   * Check if an option exists at all. Typically used for simple switches. 
   *
   * \param option String containing the sought option.
   */
  bool cmdOptionExists( const std::string &option ) const;

  /**
   * Returns the complete (reconstructed) command line.
   * @return String containing executable name and command line tokens.
   */
  std::string getCommandLine() const;

  /**
   * Returs true if all required parameters were set in the command line, false otherwise. 
   */
  bool getAllReqreuiredParametersSet() const { return this->allRequiredParametersSet; };


private:
  std::vector <std::string> tokens; ///< vector with individual command line entries
  std::string executableName;       ///< the executable name contained on the command line
  bool allRequiredParametersSet;    ///< defines if all required parameters were found
};

