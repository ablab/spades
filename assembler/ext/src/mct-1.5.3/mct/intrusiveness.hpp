// This file is part of Miscellaneous Container Templates.
//
//             https://launchpad.net/libmct/
//
// Copyright (c) 2010, 2011, 2012 Paul Pogonyshev.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.


#ifndef MCT_INTRUSIVENESS_HPP
#define MCT_INTRUSIVENESS_HPP


#include <utility>

#include <mct/impl/utils.hpp>


namespace mct
{

  // Templated structure specifying how given type could be used by its container (in
  // broad sense), i.e. "externally".  Default is no external use; it can be enabled by
  // specializing the structure.
  template <typename Type, typename = void>
  struct external_use
  { };


  template <typename Type, typename = Type>
  struct supports_external_use : impl::false_type
  { };

  template <typename Type>
  struct supports_external_use <Type, typename external_use <Type>::type> : impl::true_type
  { };


  template <typename Type, typename Value, Value Type::* field,
            bool direct_access = (impl::is_integral <Value>::value
                                  && !impl::is_same <Value, bool>::value),
            bool recurse       = supports_external_use <Value>::value>
  struct extern_use_field;

  template <typename Type, typename Value, Value Type::* field>
  struct extern_use_field <Type, Value, field, true, false>
  {
    typedef  Type   type;
    typedef  Value  value_type;

    static  const value_type&
    get (const type& structure)
    {  return structure.*field;  }

    static  void
    set (type& structure, const value_type& value)
    {  structure.*field = value;  }
  };

  template <typename Type, typename Value, Value Type::* field>
  struct extern_use_field <Type, Value, field, false, true>
  {
    typedef  Type                                       type;
    typedef  typename external_use <Value>::value_type  value_type;

    static  const value_type&
    get (const type& structure)
    {
      return external_use <Value>::get (structure.*field);
    }

    static  void
    set (type& structure, const value_type& value)
    {
      external_use <Value>::set (structure.*field, value);
    }
  };


  template <typename First, typename Second>
  struct external_use <std::pair <First, Second>,
                       typename impl::enable_if <supports_external_use <First>::value>::type>
    : extern_use_field <std::pair <First, Second>, First, &std::pair <First, Second>::first>
  { };

  template <typename First, typename Second>
  struct external_use <std::pair <First, Second>,
                       typename impl::enable_if <supports_external_use <Second>::value
                                                 && !supports_external_use <First>::value>::type>
    : extern_use_field <std::pair <First, Second>, Second, &std::pair <First, Second>::second>
  { };


  namespace impl
  {

    // By default the structure is empty.  'intrusive_storage' below makes sure to never
    // use it when it's empty, i.e. of size 1.
    template <typename type, typename = void>
    struct extern_use_wrapper
    { };


# if 0

    // Ideally we'd want this, but see comment in the preprocessor-enabled branch.
    template <typename type>
    struct extern_use_wrapper <type, typename impl::enable_if <std::is_class <type>::value>::type>
      : type
    // ...

# else

    // We currently specialize 'extern_use_wrapper' very conservatively.  At least GCC up
    // to 4.6 doesn't provide 'std::is_class', so there is no reliable way to determine if
    // 'type' is subclassable without requiring a recent compiler.  The most important
    // case, especially for maps, is 'std::pair', so we limit ourselves to that for now.
    template <typename First, typename Second>
    struct extern_use_wrapper <std::pair <First, Second> > : std::pair <First, Second>
    {
      typedef  std::pair <First, Second>  base_type;

      char  _unused;

      extern_use_wrapper&
      operator= (const extern_use_wrapper& that)
      {
        return static_cast <extern_use_wrapper&> (base_type::operator= (that));
      }

#   if MCT_CXX0X_SUPPORTED

      extern_use_wrapper&
      operator= (extern_use_wrapper&& that)
      {
        return static_cast <extern_use_wrapper&> (base_type::operator=
                                                  (std::forward <base_type> (that)));
      }

#   endif
    };

# endif

  }


  template <typename type>
  struct external_use <impl::extern_use_wrapper <type> >
    : extern_use_field <impl::extern_use_wrapper <type>, char,
                        &impl::extern_use_wrapper <type>::extern_use_wrapper::_unused>
  { };


  template <typename Value, typename = void>
  struct intrusive_storage
  {
    typedef  Value  type;
  };

  // The whole point of this type: use a wrapper that allows for external use but only if
  // 'Value' type doesn't already support external use and wrapper size is the same.
  // Also, we don't enable this specialization if size of 'extern_use_wrapper' is 1, which
  // means it was not specialized for 'Value'.
  template <typename Value>
  struct intrusive_storage <Value,
                            typename impl::enable_if <!supports_external_use <Value>::value
                                                      && sizeof (impl::extern_use_wrapper <Value>) != 1
                                                      && (sizeof (impl::extern_use_wrapper <Value>)
                                                          == sizeof (Value))>::type>
  {
    typedef  impl::extern_use_wrapper <Value>  type;
  };


  // This functionality is part of 'extern_use_field' now, so 'propagate_external_use' is
  // no longer needed and is deprecated.
  template <typename To, typename From, From To::* field>
  struct MCT_DEPRECATED propagate_external_use : extern_use_field <To, From, field>
  { };


  namespace impl
  {

    // The following wrappers over public interface are currently package-private and
    // undocumented as they don't look generally useful.  The whole purpose here is to
    // pass through constness of the first field in 'pair <const First, Second>' --- the
    // type that's stored in maps.
    //
    // Apparently, there is no way to store 'pair <K, M>' internally, but return a
    // reference to it as 'pair <const K, M>&' as required by interface (at least without
    // 'reinterpret_cast', but I'm afraid of type-punnig bugs).  So, maps do store the
    // first-field-is-const pairs, but 'const_cast' constness away when needed.  Since we
    // only need one level of this ugliness, 'hackish_external_use' is not recursive.

    template <typename Type, typename = void>
    struct hackish_external_use : external_use <Type>
    { };


    template <typename Type, typename = Type>
    struct supports_hackish_external_use : impl::false_type
    { };

    template <typename Type>
    struct supports_hackish_external_use <Type, typename hackish_external_use <Type>::type>
      : impl::true_type
    { };


    template <typename First, typename Second>
    struct hackish_external_use
      <std::pair <const First, Second>,
       typename impl::enable_if <supports_external_use <First>::value
                                 && !supports_external_use <std::pair <const First,
                                                                       Second> >::value>::type>
    {
      typedef  std::pair <const First, Second>            type;
      typedef  typename external_use <First>::value_type  value_type;

      static  const value_type&
      get (const type& pair)
      {
        return external_use <First>::get (pair.first);
      }

      static  void
      set (type& pair, const value_type& value)
      {
        external_use <First>::set (const_cast <First&> (pair.first), value);
      }
    };

  }

}


#endif  // Multi-inclusion guard.


// Local variables:
// mode: c++
// c-basic-offset: 2
// indent-tabs-mode: nil
// fill-column: 90
// End:
