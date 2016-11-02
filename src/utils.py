"""
Some useful tools.

Author:       Benedikt J. Daurer (benedikt@xray.bmc.uu.se)
Last update:  March 5, 2014
"""
import os, sys, datetime, time
import numpy as np
__all__ = ['progressbar']

# ----------------------------
# GET WIDTH/HEIGHT OF TERMINAL
# ----------------------------
def tty_size():
    """
    Get the width and height of the terminal in number of characters

    Parameters
    ----------
    default : integer, optional
        Value that will be returned if the function fails; the default is 80.

    Returns
    -------
    tty_size : tuple
        The width and height of the tty.

    Reference
    ---------
    http://stackoverflow.com/questions/566746/how-to-get-console-window-width-in-python
    """
    env = os.environ
    def ioctl_GWINSZ(fd):
        try:
            import fcntl, termios, struct, os
            cr = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ,
        '1234'))
        except:
            return
        return cr
    cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if not cr:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            cr = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass
    if not cr:
        cr = (env.get('LINES', 25), env.get('COLUMNS', 80))
    return int(cr[1]), int(cr[0])
    

# -------------------------
# PRINT A TEMPORARY MESSAGE
# -------------------------
def print_tmp(message = '', align='l'):
    """
    Print a message that will be overwritten by any subsequent output.

    Parameters
    ----------
    message : string, optional
        The message to be printed.
        If the string is empty (default) it clears the line.
        align : character, optional
        The alignment of the message; left by default.
        This value can be one of 'l' (left), 'c' (center) or 'r' (right).

    Returns
    -------    
    tmp_message : string
        The printed message extended to the with of the tty.

    See Also
    --------
    tty_size

    Example
    -------
    >>> import time
    >>> print_tmp('Hello World!')
    >>> print_tmp('>>> Short message <<<', 'c'); time.sleep(1.); print_tmp()
    >>> print_tmp('Good Bye!', 'r')
    """
    
    # only take first character of the align parameter
    align = align[0]
    
    # get width of the terminal
    tty_width = tty_size()[0]
    
    message_len = len(message)
    if message_len >= tty_width:
        # limit message length to terminal width
        message = message[:tty_width]
    elif align == 'c':
        message = message.rjust((tty_width + message_len) / 2).ljust(tty_width)
    elif align == 'r':
        message = message.rjust(tty_width)
    else:  # default, if unknown alignment
        message = message.ljust(tty_width)
    
    sys.stdout.write("\r" + message)
    if message == '':
        sys.stdout.write("\r")  # clear the line
    sys.stdout.flush()
    
    return message


def progressbar(progress = 1., name = '', start_timestamp = None, bar_char = '=', tty_width = None):
    """
    Print a progressbar.

    Parameters
    ----------
    progress : float, optional
        The progress as value between 0 and 1.
        If progress is 1. (default), a linebreak will be appended to the output.
    name : string, optional
        A name displayed for the progress; empty by default.
    start_timestamp : timestamp, optional
        Estimate and display the remaining time from the time given.
        This feature off by default.
    bar_char : character, optional
        The character to build up the progressbar, '=' by default.
    tty_width : integer, optional
        The total number of characters for the progressbar. This defaults to the
        width of the tty.

    Returns
    -------    
    tmp_message : string
        The printed progressbar.

    See Also
    --------
    tty_size
    print_tmp

    Example
    -------
    >>> import utils
    >>> import time
    >>> N = 10
    >>> start = time.time()
    >>> for i in range(N):
    ...     progress = float(i) / N
    ...     utils.progressbar(progress, 'do something', start)
    ...     utils.pause(.5)
    >>> utils.progressbar(1., 'do something', start)
    """
    # enforce progressbar to be in [1,0] and bar_char to be only one character
    progress = np.clip(progress, a_min=0., a_max=1.)
    bar_char = bar_char[0]
    
    # get width of the terminal
    if tty_width is None:
        tty_width = tty_size()[0]
    
    head = name + ': [' if name != '' else '['               # name
    tail = '] ' + str(int(progress * 100.)).rjust(3) + '% '  # percent
    
    if start_timestamp is not None:                          # remaining time
        if progress != 0.:
            if progress != 1.:
                tail += '- ' + str(datetime.timedelta(0, int((time.time() - start_timestamp) * (1. / progress - 1.)))) + ' '
            else:                                            # total time here
                tail += '- ' + str(datetime.timedelta(0, int(time.time() - start_timestamp))) + ' '
        else:
            tail += '- ?:??:?? '

    bar_len = tty_width - len(head) - len(tail)
    if progress != 1.:
        if progress != 0.:
            output = head + \
                (bar_char * tty_width + '>')[-int(np.ceil(progress * bar_len)):].ljust(bar_len) + \
                tail
        else:
            output = head + ' ' * bar_len + tail
        sys.stdout.write("\r" + output)
    else:
        output = head + \
            ("==== finished ".ljust(bar_len, '=') if bar_len > 14 else '=' * bar_len) + \
            tail
        sys.stdout.write("\r" + output + os.linesep)
    sys.stdout.flush()
    
    return output
