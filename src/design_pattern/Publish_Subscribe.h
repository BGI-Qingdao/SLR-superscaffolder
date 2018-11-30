#ifndef __DESIGN_PATTERN_PUBLISH_SUBSCRIBE_H__
#define __DESIGN_PATTERN_PUBLISH_SUBSCRIBE_H__

#include <set>
#include <cassert>

namespace BGIQD {
    namespace DesignPattern {

        template<class Message>
            struct ISubscriber
            {
                virtual void update_msg(const Message & item) = 0;
            };

        template<class M, class S = ISubscriber<M> >
            struct IPublisher
            {
                typedef M Message ;

                typedef S Subscriber ;

                void add( Subscriber * ptr )
                {
                    assert(ptr != NULL );
                    _subsribers.insert(ptr);
                }

                void remove ( Subscriber * ptr )
                {
                    _subsribers.erase(ptr);
                }

                void set_msg( const Message & msg )
                {
                    _msg = msg ;
                }

                void notify() 
                {
                    for( auto x : _subsribers )
                    {
                        assert( x != NULL );
                        x->update_msg(_msg);
                    }
                };

                void notify_with_msg( const Message & msg )
                {
                    for( auto x : _subsribers )
                    {
                        assert( x != NULL );
                        x->update_msg(msg);
                    }
                }

                Message get_msg() const 
                {
                    return _msg ;
                }

                Message _msg;
                std::set<Subscriber*> _subsribers;
            };
    }
}

#endif //__DESIGN_PATTERN_PUBLISH_SUBSCRIBE_H__
