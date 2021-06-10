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

/** Structure holding the basic information about a single command line option. 
 */
struct CommandLineOption{
  unsigned int numberOfComponents; ///< Number of elements expected after this flag on the command line. A switch will be 0, a single item 1, etc. 
  bool required;                   ///< Indicating if this parameters is essential to run the program. If a given parameter marked by required is not found on the command line a warning will be printed. 
  std::string description;         ///< Detailed description of the parameter. This will be used to generate the help message on the command line. 
  std::string format;              ///< String describing the format of the parameter, e.g. <int>, <filename>, etc.
};


/** Function to split a string by a delimiter into a vector of strings. Items will only be added if they are not empty. 
 *  I.e. if the delimiter is a comma a string containing two commas will not add an empty element. 
 * 
 *  \param stringToSplitIn The string which will be split by the given delimiter
 *  \param delimiter       Delimiter which will be used to split the string
 */
std::vector<std::string> splitStringbyDelimiter(const std::string & stringToSplitIn, const std::string & delimiter);


/** Function to print a selection of parameters. 
 *  The parameters are assumed to belong to a semantic group of items such as input and put. This can be specified in the 
 *  section heading. The final format will be looking like this:
 *  \code
 *    Section heading
 *    ~~~~~~~~~~~~~~~
 *    -opt1 <format> The description will go here and will break after
 *                   the maximum width was reached. 
 *    -opt2 <format> Alterantively a backslash-n
 *                   can be inserted in the descriptionb to force a line
 *                   break. 
 * 
 *    |<----------------------------maxWidth---------------------------->|
 *  \endcode
 * 
 *  \param commandLineOptions The map of all commandline options. The format and description will be used to generate the output. 
 *  \param optionSectionHeading A string that will be printed as a heading above the parameters
 *  \param maxWidth The total width of the line printed. A warning only may be printed if the width is not sufficient to accomodate 
 *                  the option and the format (i.e. no space left for the description). The maximum width may be exceeded if the first 
 *                  word of the description is already too long. 
 *  \param optionsToPrint A vector of option keys to be looked up in the given commandLineOptions. If present the format and description
 *                        of the commandLineOptions will be used to generate the output as shown above. 
 */
void printFormattedCommandLineOptions(const std::map<std::string, CommandLineOption>& commandLineOptions, 
  const std::string& optionSectionHeading, const unsigned int maxWidth, const std::vector<std::string>& optionsToPrint );


/** Class implementing a simple command line parser
 *  Command line parser (heavily) adapted from 
 *  http://stackoverflow.com/questions/865668/ddg#868894
 *  credit to iain.
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
  std::vector <std::string> tokens; ///< Vector with individual command line entries
  std::string executableName;       ///< The executable name contained on the command line
  bool allRequiredParametersSet;    ///< Indicates if all required parameters were found. Use the getter funciton \ref getAllReqreuiredParametersSet().
};

