// -*- C++ -*-

#ifndef CONF_MAN_HH
#define CONF_MAN_HH

#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

#include <globals.hh>

//_____________________________________________________________________________
class ConfMan
{
public:
  static G4String  ClassName( void );
  static ConfMan& GetInstance( void );
  ~ConfMan( void );

private:
  ConfMan( void );
  ConfMan( const ConfMan& );
  ConfMan& operator=( const ConfMan& );

private:
  typedef std::map<G4String, G4String> StrList;
  typedef std::map<G4String, G4double> DoubleList;
  typedef std::map<G4String, G4int>    IntList;
  typedef std::map<G4String, G4bool>   BoolList;
  G4String     m_conf_key;
  G4String     m_conf_dir;
  G4bool       m_is_ready;
  StrList      m_file;
  StrList      m_string;
  DoubleList   m_double;
  IntList      m_int;
  BoolList     m_bool;

public:
  // G4bool    Finalize( void );
  // G4bool    FinalizeProcess( void );
  template <typename T>
  static const T& Get( const G4String& key );
  G4bool    Initialize( void );
  G4bool    Initialize( const G4String& file_name );
  G4bool    Initialize( const G4String& file_name ,const G4String& file_name2 );
  G4bool    InitializeHistograms( void );
  G4bool    InitializeParameterFiles( void );
  template <typename T>
  G4bool    InitializeParameter( void );
  template <typename T>
  G4bool    InitializeParameter( const G4String& key );
  template <typename T>
  G4bool    InitializeParameter( const G4String& key1,
				 const G4String& key2 );
  G4bool    IsReady( void ) const { return m_is_ready; }

private:
  G4String FilePath( const G4String& src ) const;
  G4bool  ShowResult( G4bool status, const G4String& name ) const;
};

//_____________________________________________________________________________
inline G4String
ConfMan::ClassName( void )
{
  static const G4String s_name("ConfMan");
  return s_name;
}

//_____________________________________________________________________________
inline ConfMan&
ConfMan::GetInstance( void )
{
  static ConfMan s_instance;
  return s_instance;
}

//_____________________________________________________________________________
template <>
inline const G4String&
ConfMan::Get<G4String>( const G4String& key )
{
  return GetInstance().m_string[key];
}

//_____________________________________________________________________________
template <>
inline const G4double&
ConfMan::Get<G4double>( const G4String& key )
{
  return GetInstance().m_double[key];
}

//_____________________________________________________________________________
template <>
inline const G4int&
ConfMan::Get<G4int>( const G4String& key )
{
  return GetInstance().m_int[key];
}

//_____________________________________________________________________________
template <>
inline const G4bool&
ConfMan::Get<G4bool>( const G4String& key )
{
  return GetInstance().m_bool[key];
}

//_____________________________________________________________________________
inline G4bool
ConfMan::ShowResult( G4bool status, const G4String& name ) const
{
  if( status )
    std::cout << std::setw(20) << std::left
	      << " ["+name+"]"
	      << "-> Initialized" << std::endl;
  else
    std::cout << std::setw(20) << std::left
	      << " ["+name+"]"
	      << "-> Failed" << std::endl;
  return status;
}

//_____________________________________________________________________________
template <typename T>
inline G4bool
ConfMan::InitializeParameter( void )
{
  return
    ShowResult( T::GetInstance().Initialize(),
		T::GetInstance().ClassName() );
}

//_____________________________________________________________________________
template <typename T>
inline G4bool
ConfMan::InitializeParameter( const G4String& key )
{
  return
    ShowResult( T::GetInstance().Initialize(m_file[key]),
		T::GetInstance().ClassName() );
}

//_____________________________________________________________________________
template <typename T>
inline G4bool
ConfMan::InitializeParameter( const G4String& key1,
			      const G4String& key2 )
{
  return
    ShowResult( T::GetInstance().Initialize(m_file[key1],
					    m_file[key2]),
		T::GetInstance().ClassName() );
}

#endif
